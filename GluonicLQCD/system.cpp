#include <random>   // For Mersenne-Twister19937
#include <ctime>    // For random seed
#include <cmath>    // For exp()
#include <fstream>
#include <iostream>
#include <cstdio>   // For io C-style handling.
#include <cstdlib>
#include <mpi.h>
#include "system.h"
#include "actions/action.h"
#include "correlators/correlator.h"
#include "functions.h"
#include "links.h"
#include "matrices/su3matrixgenerator.h"
#include "parallelization/neighbourlist.h"
#include "parallelization/neighbours.h"
#include "parallelization/indexorganiser.h"


using std::cout;
using std::endl;

System::System(int NSpatial, int NTemporal, int NCf, int NCor, int NTherm, double beta, double seed, Correlator *correlator, Action *S, int numprocs, int processRank)
{
    /*
     * Class for calculating correlators using the System algorithm.
     * Takes an action object as well as a Gamma functional to be used in the action.
     */
    m_NSpatial = NSpatial; // Spatial dimensions
    m_NTemporal = NTemporal; // Time dimensions
    m_latticeSize = NSpatial*NSpatial*NSpatial*NTemporal;
    m_NCf = NCf; // Number of configurations to run for
    m_NCor = NCor;
    m_NTherm = NTherm;
    m_beta = beta;
    m_numprocs = numprocs;
    m_processRank = processRank;
    setAction(S);
    setCorrelator(correlator);
    m_GammaPreThermalization = new double[m_NTherm*m_NCor+1];
    m_Gamma = new double[m_NCf]; // Correlator values
    m_GammaSquared = new double[m_NCf];

    // For parallelization
    m_neighbourLists = new Neighbours;
    m_indexHandler = new IndexOrganiser(m_processRank);

    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uni_dist(0,1);
    m_generator = gen;
    m_uniform_distribution = uni_dist;
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        m_Gamma[alpha] = 0;
        m_GammaSquared[alpha] = 0;
    }
}

System::~System()
{
    /*
     * Class destructor
     */
    delete [] m_lattice;
    delete [] m_Gamma;
    delete [] m_GammaSquared;
    delete [] m_GammaPreThermalization;
}

void System::subLatticeSetup()
{
    /*
     * Sets up the sub-lattices. Adds +2 in every direction to account for sharing of s.
     */
    if (m_numprocs % 2 != 0) {
        cout << "Error: odd number of processors --> exiting." << endl;
        exit(1);
    }
    int restProc = m_numprocs;

    // Sets up sub lattice dimensionality without any splitting
    for (int i = 0; i < 3; i++) {
        m_NTrue[i] = m_NSpatial;
    }
    m_NTrue[3] = m_NTemporal;

    // Iteratively finds and sets the sub-lattice cube sizes
    while (restProc >= 2) {
        for (int i = 0; i < 4; i++) { // Conts from x to t
            m_NTrue[i] /= 2;
            restProc /= 2;
            if (restProc < 2) break;
        }
    }
    m_subLatticeSize = 1;
    for (int i = 0; i < 4; i++) {
        m_subLatticeSize *= m_NTrue[i]; // Gets the total size of the sub-lattice(without faces)
    }
    m_lattice = new Links[m_subLatticeSize];

    // Sets up number of processors per dimension
    for (int i = 0; i < 3; i++) {
        m_processorsPerDimension[i] = m_NSpatial / m_NTrue[i];
    }
    m_processorsPerDimension[3] = m_NTemporal / m_NTrue[3];
    m_neighbourLists->initialize(m_processRank, m_numprocs, m_processorsPerDimension);

    m_indexHandler->setN(m_NTrue);
    m_indexHandler->setNeighbourList(m_neighbourLists);
    m_S->initializeIndexHandler(m_indexHandler);
    m_correlator->initializeIndexHandler(m_indexHandler);
    m_S->setN(m_NTrue);
    m_correlator->setN(m_NTrue);
    m_correlator->setLatticeSize(m_subLatticeSize);
}

void System::latticeSetup(SU3MatrixGenerator *SU3Generator, bool hotStart)
{
    /*
     * Sets up the lattice and its matrices.
     */
    subLatticeSetup();

    m_SU3Generator = SU3Generator;
    if (hotStart) {
        // All starts with a completely random matrix.
        for (int i = 0; i < m_latticeSize; i++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                m_lattice[i].U[mu] = m_SU3Generator->generateRandom(); // Fully random
//            m_lattice[i].U[mu] = m_SU3Generator->generateRST(); // Random close to unity
            }
        }
    } else {
        // Cold start: everything starts out at unity.
        for (int x = 0; x < m_NTrue[0]; x++) {
            for (int y = 0; y < m_NTrue[1]; y++) {
                for (int z = 0; z < m_NTrue[2]; z++) {
                    for (int t = 0; t < m_NTrue[3]; t++) {
                        for (int mu = 0; mu < 4; mu++) {
                            m_lattice[getIndex(x,y,z,t,m_NTrue[1],m_NTrue[2],m_NTrue[3])].U[mu] = m_SU3Generator->generateIdentity();
                        }
                    }
                }
            }
        }
    }
    if (m_processRank == 0) {
        cout << "Lattice setup complete" << endl;
    }
}

void System::updateLink(int latticeIndex, int mu)
{
    /*
     * Private function used for updating our system. Updates a single gauge link.
     * Arguments:
     *  i   : spacetime getIndex
     *  mu  : Lorentz getIndex
     */
//    SU3 X = m_SU3Generator->generateRandom(); // Generates a random matrix, SHOULD BE MODIFIED TO X = RST, page 83 Gattinger & Lang
//    SU3 X = m_SU3Generator->generateRST();
//    m_updatedMatrix = X*m_lattice[latticeIndex].U[mu];
    m_updatedMatrix = m_SU3Generator->generateRST()*m_lattice[latticeIndex].U[mu]; // Shorter method of updating matrix
}

void System::update()
{
    /*
     * Sweeps the entire Lattice, and gives every matrix a chance to update.
     */
    for (int x = 0; x < m_NTrue[0]; x++) {
        for (int y = 0; y < m_NTrue[1]; y++) {
            for (int z = 0; z < m_NTrue[2]; z++) {
                for (int t = 0; t < m_NTrue[3]; t++) {
                    for (int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice, x, y, z, t, mu);
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(getIndex(x, y, z, t, m_NTrue[1], m_NTrue[2], m_NTrue[3]), mu);
                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu);
                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[getIndex(x, y, z, t, m_NTrue[1], m_NTrue[2], m_NTrue[3])].U[mu].copy(m_updatedMatrix);
                                m_acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
}


void System::runMetropolis(bool storePreObservables, bool storeConfigs)
{
//    loadFieldConfiguration("conf0.bin");
    // Variables for checking performance of the update.
    clock_t preUpdate, postUpdate;
    double updateStorer = 0;

    // Calculating correlator before any updates have began.
    m_GammaPreThermalization[0] = m_correlator->calculate(m_lattice);

    // Summing and sharing correlator to all processors before any updates has begun
    MPI_Allreduce(&m_GammaPreThermalization[0], &m_GammaPreThermalization[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Dividing by the number of processors in order to get the correlator.
    m_GammaPreThermalization[0] /= double(m_numprocs);
    if (m_processRank == 0) cout << "Pre-thermialization correlator:  " << m_GammaPreThermalization[0] << endl;

    // Running thermalization
    for (int i = 0; i < m_NTherm*m_NCor; i++)
    {
        preUpdate = clock();
        update();
        postUpdate = clock();
        updateStorer += ((postUpdate - preUpdate)/((double)CLOCKS_PER_SEC));
        if ((i-1) % 20 == 0) {
            if (m_processRank == 0) {
                cout << "Avg. time per update every 20th update: " << updateStorer/(i+1) << " sec" << endl;
            }
        }
        postUpdate = clock();
        // ========= REDO THIS =========
        // Print correlator every somehting or store them all(useful when doing the thermalization).
        if (storePreObservables) {
            // Calculating the correlator
            m_GammaPreThermalization[i+1] = m_correlator->calculate(m_lattice);
            // Summing and sharing results across the processors
            MPI_Allreduce(&m_GammaPreThermalization[i+1], &m_GammaPreThermalization[i+1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // Averaging the results
            m_GammaPreThermalization[i+1] /= double(m_numprocs);
            if (m_processRank == 0) cout << i << " " << m_GammaPreThermalization[i+1] << endl; // Printing evolution of system
        } else if ((i+1) % 10 == 0) {
            m_GammaPreThermalization[int((i+1)/10.0)+1] = m_correlator->calculate(m_lattice);
            MPI_Allreduce(&m_GammaPreThermalization[int((i+1)/10.0)+1], &m_GammaPreThermalization[int((i+1)/10.0)+1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            m_GammaPreThermalization[int((i+1)/10.0)+1] /= double(m_numprocs);
        }
    }

    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    m_acceptanceCounter = int(double(m_acceptanceCounter)/double(m_numprocs));

    // Printing post-thermalization correlator and acceptance rate
    if (m_processRank == 0) {
        cout << "Post-thermialization correlator: " << m_GammaPreThermalization[m_NTherm*m_NCor] << endl;
        cout << "Termalization complete. Acceptance rate: " << m_acceptanceCounter/double(4*m_latticeSize*m_nUpdates*m_NTherm*m_NCor) << endl;
    }

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;

    // Main part of algorithm
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
        }
        // Averaging the gamma values
        m_Gamma[alpha] = m_correlator->calculate(m_lattice);
        MPI_Allreduce(&m_Gamma[alpha], &m_Gamma[alpha], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        m_Gamma[alpha] /= double(m_numprocs);
        // Writing to file
        if (m_processRank == 0) writeConfigurationToFile();
    }
    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    m_acceptanceCounter = int(double(m_acceptanceCounter)/double(m_numprocs));
    if (m_processRank == 0) cout << "System completed." << endl;
}

void System::sampleSystem()
{
    /*
     * For sampling statistics/getting correlators
     */
    cout << "Not implemented yet." << endl;
}

void System::runBasicStatistics()
{
    /*
     * Class instance for sampling statistics from our system.
     */
    double averagedGammaSquared = 0;
    // Performing an average over the Monte Carlo obtained values
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        m_averagedGamma += m_Gamma[alpha];
        averagedGammaSquared += m_Gamma[alpha]*m_Gamma[alpha];
    }
    m_averagedGamma /= double(m_NCf);
    averagedGammaSquared /= double(m_NCf);
    m_varianceGamma = (averagedGammaSquared - m_averagedGamma*m_averagedGamma)/double(m_NCf);
    m_stdGamma = sqrt(m_varianceGamma);
    if (m_processRank == 0) cout << m_averagedGamma << ", std = " << m_stdGamma << ", variance = " << m_varianceGamma << endl;
}

void System::writeDataToFile(std::string filename, bool preThermalizationGamma)
{
    /*
     * For writing the raw Gamma data to file.
     * Arguments:
     *  filename                : to write stats to
     *  preThermalizationGamma  : if we are writing the gama
     */
    if (m_processRank == 0) {
        for (int i = 0; i < m_NCf; i++) {
            cout << m_Gamma[i] << endl;
        }
        std::ofstream file;
        file.open(filename + ".dat");
        file << "acceptanceCounter " << getAcceptanceRate() << endl;
        file << "NCor " << m_NCor << endl;
        file << "NCf " << m_NCf << endl;
        file << "NTherm " << m_NTherm << endl;
        file << "AverageGamma " << m_averagedGamma << endl;
        file << "VarianceGamma " << m_varianceGamma << endl;
        file << "stdGamma " << m_stdGamma << endl;
        if (preThermalizationGamma) {
            for (int i = 0; i < m_NTherm*m_NCor+1; i++) {
                file << m_GammaPreThermalization[i] << endl;
            }
        }
        for (int i = 0; i < m_NCf; i++) {
            file << m_Gamma[i] << endl;
        }
        file.close();
        cout << filename << " written" << endl;
    }
}


void System::printAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the System algorithm.
     */
    if (m_processRank == 0) printf("Acceptancerate: %.16f \n", getAcceptanceRate());
}

double System::getAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the System algorithm.
     */
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_nUpdates*m_latticeSize*4); // Times 4 from the Lorentz indices
}

void System::writeConfigurationToFile()
{
    /*
     * C-method for writing out configuration to file.
     * Arguments:
     * - filename
     */
    // IMPLEMENT FULL WRITE TO FILE FOR LATTICE HERE
//    FILE *file; // C method
//    file = fopen((m_outputFolder + "/" + m_filename + "_beta" + std::to_string(m_beta) + "_config" + std::to_string(configNumber) + ".bin").c_str(), "wb");

    // buffer?

    MPI_Datatype view;
    // TRYING TO IMPLEMENT MPI_CREATE_SUB_LATTICE

    for (int t = 0; t < m_NTemporal; t++) {
        for (int z = 0; z < m_NSpatial; z++) {
            for (int y = 0; y < m_NSpatial; y++) {
                for (int x = 0; x < m_NSpatial; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fwrite(&m_lattice[getIndex(x, y, z, t, m_NTrue[1], m_NTrue[2], m_NTrue[3])].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }

//    // OLD WRITE TO FILE(SCALAR)
//    for (int t = 0; t < m_NTrue[3]; t++) {
//        for (int z = 0; z < m_NTrue[2]; z++) {
//            for (int y = 0; y < m_NTrue[1]; y++) {
//                for (int x = 0; x < m_NTrue[0]; x++) {
//                    for (int mu = 0; mu < 4; mu++) {
//                        fwrite(&m_lattice[getIndex(x, y, z, t, m_NTrue[1], m_NTrue[2], m_NTrue[3])].U[mu],sizeof(SU3),1,file);
//                    }
//                }
//            }
//        }
//    }
//    fclose(file);
    cout << m_outputFolder + m_outputFolder + "/" + m_filename + "_beta" + std::to_string(m_beta) + "_config" + std::to_string(configNumber) + ".bin"  + " written" << endl;
}

void System::loadFieldConfiguration(std::string filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */
    // IMPLEMENT LOAD FROM FILE FOR LATTICE HERE
    FILE *file; // C method
    file = fopen((m_inputFolder +"_p" + std::to_string(m_processRank) + filename).c_str(), "rb");
    for (int t = 0; t < m_NTrue[3]; t++) {
        for (int z = 0; z < m_NTrue[2]; z++) {
            for (int y = 0; y < m_NTrue[1]; y++) {
                for (int x = 0; x < m_NTrue[0]; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fread(&m_lattice[getIndex(x, y, z, t, m_NTrue[1], m_NTrue[2], m_NTrue[3])].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }
    fclose(file);
    cout << m_inputFolder + filename  + " loaded" << endl;
}
