#include <random>   // For Mersenne-Twister19937
#include <ctime>    // For random seed
#include <cmath>    // For exp()
#include <fstream>
#include <iostream>
#include <cstdio>   // For io C-style handling.
#include <cstdlib>
#include "system.h"
#include "actions/action.h"
#include "correlators/correlator.h"
#include "functions.h"
#include "links.h"
#include "matrices/su3matrixgenerator.h"
#include "mpi.h"

//TEMP
#include "unittests.h"

using std::cout;
using std::endl;

System::System(int N, int N_T, int NCf, int NCor, int NTherm, double a, double L, double seed, Correlator *correlator, Action *S, int numprocs, int processRank)
{
    /*
     * Class for calculating correlators using the System algorithm.
     * Takes an action object as well as a Gamma functional to be used in the action.
     */
    m_N = N; // Spatial dimensions
    m_N_T = N_T; // Time dimensions
    m_latticeSize = N*N*N*N_T;
    m_NCf = NCf; // Number of configurations to run for
    m_NCor = NCor;
    m_NTherm = NTherm;
    m_a = a;
    m_L = L;
    m_numprocs = numprocs;
    m_processRank = processRank;
    setAction(S);
    setCorrelator(correlator);
//    m_lattice = new Links[m_latticeSize]; // Lattice, contigious memory allocation
//    m_GammaPreThermalization = new double[m_NTherm*m_NCor/10];
    m_GammaPreThermalization = new double[m_NTherm*m_NCor+1];
    m_Gamma = new double[m_NCf]; // Correlator values
    m_GammaSquared = new double[m_NCf];

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
}

void System::subLatticeDimensionsSetup()
{
    /*
     * Sets up the sub-lattices. Adds +2 in every direction to account for sharing of phases.
     */
    if (m_numprocs % 2 != 0) {
        cout << "Error: odd number of processors --> exiting." << endl;
        exit(1);
    }
    int restProc = m_numprocs;

    // Sets up sub lattice dimensionality without any splitting
    for (int i = 0; i < 3; i++) {
        m_subLatticeDimensions[i] = m_N;
    }
    m_subLatticeDimensions[3] = m_N_T;

//    // TEST==========================================================
//    if (m_processRank == 0) {
//        cout << "Processor: " << m_processRank << endl;
//        for (int i = 0; i < 4; i++) {
//            cout << m_subLatticeDimensions[i] << endl;
//        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    // ==============================================================

    // Iteratively finds the lattice cube sizes
    while (restProc >= 2) {
//        for (int i = 3; i > 0; i--) {
        for (int i = 0; i < 4; i++) {
            m_subLatticeDimensions[i] /= 2;
            restProc /= 2;
            if (restProc < 2) break;
        }
    }
    m_subLatticeSize = 1;
    m_trueSubLatticeSize = 1;
    for (int i = 0; i < 4; i++) {
        m_subLatticeSize *= m_subLatticeDimensions[i]; // Gets the total size of the sub-lattice(without phases)
        m_trueSubLatticeDimensions[i] = m_subLatticeDimensions[i] + 2; // Adds a phase
        m_trueSubLatticeSize *= m_trueSubLatticeDimensions[i]; // Gets the total size of the sub-lattice(with phases)
    }
    m_latticeSize = m_subLatticeSize; // OK TO DO THIS? CREATE GLOBAL SUB LATTICE CONSTANTS ECT?
    m_lattice = new Links[m_trueSubLatticeSize];

    cout << "Process rank: " << m_processRank << endl;
    cout << "m_subLatticeSize = " << m_subLatticeSize << endl;
    cout << "m_trueSubLatticeSize = " << m_trueSubLatticeSize << endl;

    // Sets up number of processors per dimension
    for (int i = 0; i < 3; i++) {
        m_processorsPerDimension[i] = m_N / m_subLatticeDimensions[i];
    }
    m_processorsPerDimension[3] = m_N_T / m_subLatticeDimensions[3];

    // TODO: Create function that takes Npx,Npy,Npz,Npt as arguments(that is, processors in each direction)

    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */

    // Sets up neighbour lists
//    m_neighbourList = new int[8];
//    m_allNeighbourLists = new int[8*m_numprocs];
//    neighbourIndex(m_processRank, , m_numprocs);
//    fillNeighbourLists();

    // Assign x direction neighbours
    // Assign y direction neighbours
    // Assign z direction neighbours
    // Assign t direction neighbours

    // PRINTS ===========================================================
//    cout << "Processor: " << m_processRank << endl;
//    for (int i = 0; i < 4; i++) {
//        cout << "m_dim     = " << m_subLatticeDimensions[i] << endl;
//        cout << "m_trueDim = " << m_trueSubLatticeDimensions[i] << endl;
//    }
//    // ==================================================================

//    cout << "EXITS AT SUB LATTICE DIM SETTUP" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    exit(1);
}

void System::latticeSetup(SU3MatrixGenerator *SU3Generator, bool hotStart)
{
    /*
     * Sets up the lattice and its matrices.
     */
    // PARALLELIZE HERE?? TIME IT!
    subLatticeDimensionsSetup();
    // Also, set up blocks to use!!
    m_SU3Generator = SU3Generator;
    if (hotStart) {
        // All starts with a completely random matrix.
//        for (int i = 0; i < m_latticeSize; i++)
        for (int i = 0; i < m_subLatticeSize; i++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                m_lattice[i].U[mu] = m_SU3Generator->generateRandom();
//            m_lattice[i].U[mu] = m_SU3Generator->generateRST();
            }
        }
    } else {
//        for (int i = 0; i < m_latticeSize; i++)
        for (int i = 0; i < m_subLatticeSize; i++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                // All starts with a completely random matrix. Observable should be 1.
                m_lattice[i].U[mu] = m_SU3Generator->generateIdentity();
            }
        }
    }
}

void System::updateLink(int latticeIndex, int mu)
{
    /*
     * Private function used for updating our system. Updates a single gauge link.
     * Arguments:
     *  i   : spacetime index
     *  mu  : Lorentz index
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
    // PARALLELIZE HERE
    for (int x = 0; x < m_N; x++) {
        for (int y = 0; y < m_N; y++) {
            for (int z = 0; z < m_N; z++) {
                for (int t = 0; t < m_N_T; t++) {
                    for (int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice, x, y, z, t, mu);
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(index(x, y, z, t, m_N, m_N_T), mu);
                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu);
                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[index(x, y, z, t, m_N, m_N_T)].U[mu].copy(m_updatedMatrix);
                            }
                            else
                            {
                                m_acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
}


void System::runMetropolis(bool storePreObservables)
{
//    loadFieldConfiguration("conf0.bin");
    m_GammaPreThermalization[0] = m_correlator->calculate(m_lattice);
    cout << "Pre-thermialization correlator:  " << m_GammaPreThermalization[0] << endl;
    // Running thermalization
    for (int i = 0; i < m_NTherm*m_NCor; i++)
    {
        update();
        // Print correlator every somehting or store them all(useful when doing the thermalization)
        if (storePreObservables) {
            m_GammaPreThermalization[i+1] = m_correlator->calculate(m_lattice);
        } else if ((i+1) % 10 == 0) {
            m_GammaPreThermalization[int((i+1)/10.0)+1] = m_correlator->calculate(m_lattice);
        }
    }
    cout << "Post-thermialization correlator: " << m_GammaPreThermalization[m_NTherm*m_NCor + 1] << endl;
    cout << "Termalization complete. Acceptance rate: " << m_acceptanceCounter/double(4*m_latticeSize*m_nUpdates*m_NTherm*m_NCor) << endl;
    delete [] m_GammaPreThermalization; // De-allocating as this is not needed anymore

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;
    // Main part of algorithm
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
        }
        m_Gamma[alpha] = m_correlator->calculate(m_lattice);
    }
    cout << "System completed." << endl;
    cout << m_correlator->calculate(m_lattice) << endl;
}

void System::sampleSystem()
{
    /*
     * For sampling statistics/getting correlators
     */
    cout << "Not implemented yet." << endl;
}

void System::getStatistics()
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
    cout << m_averagedGamma << ", std = " << m_stdGamma << ", variance = " << m_varianceGamma << endl;
}

void System::writeDataToFile(std::string filename, bool preThermalizationGamma)
{
    /*
     * For writing the raw Gamma data to file.
     * Arguments:
     *  filename                : to write stats to
     *  preThermalizationGamma  : if we are writing the gama
     */
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
//        for (int i = 0; i < m_NTherm*m_NCor/10; i++) {
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


void System::printAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the System algorithm.
     */
    printf("Acceptancerate: %.16f \n", getAcceptanceRate()); // Times 4 from the Lorentz indices
}

double System::getAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the System algorithm.
     */
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_nUpdates*m_latticeSize*4);
}

void System::writeConfigurationToFile(std::string filename)
{
    /*
     * C-method for writing out configuration to file.
     * Arguments:
     * - filename
     */
    FILE *file; // C method
    file = fopen((m_outputFolder+filename).c_str(), "wb");
    for (int t = 0; t < m_N_T; t++) {
        for (int z = 0; z < m_N; z++) {
            for (int y = 0; y < m_N; y++) {
                for (int x = 0; x < m_N; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fwrite(&m_lattice[index(x,y,z,t,m_N,m_N_T)].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }
    fclose(file);
    cout << m_outputFolder + filename  + " written" << endl;
}

void System::loadFieldConfiguration(std::string filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */
    FILE *file; // C method
    file = fopen((m_inputFolder + filename).c_str(), "rb");
    for (int t = 0; t < m_N_T; t++) {
        for (int z = 0; z < m_N; z++) {
            for (int y = 0; y < m_N; y++) {
                for (int x = 0; x < m_N; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fread(&m_lattice[index(x,y,z,t,m_N,m_N_T)].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }
    fclose(file);
    cout << m_inputFolder + filename  + " loaded" << endl;
}