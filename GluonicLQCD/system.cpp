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

System::System(int NSpatial, int NTemporal, int NCf, int NCor, int NTherm, int NUpdates, double beta, double seed, Correlator *correlator, Action *S, int numprocs, int processRank)
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
    m_NUpdates = NUpdates;
    m_beta = beta;
    m_numprocs = numprocs;
    m_processRank = processRank;
    setAction(S);
    setCorrelator(correlator);
    m_GammaPreThermalization = new double[m_NTherm+1];
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
    else if (m_numprocs % 16 != 0) {
        if (m_processRank == 0) cout << "Warning: running for an number of prosessors not divisible by 16 may slow down performance due to the geometry of the lattice." << endl;
    }
    int restProc = m_numprocs;

    // Sets up sub lattice dimensionality without any splitting
    for (int i = 0; i < 3; i++) {
        m_N[i] = m_NSpatial;
    }
    m_N[3] = m_NTemporal;

    // Iteratively finds and sets the sub-lattice dimensions
    while (restProc >= 2) {
        for (int i = 0; i < 4; i++) { // Counts from x to t
//        for (int i = 4; i >= 0; i--) { // Counts from x to t
            m_N[i] /= 2;
            restProc /= 2;
            if (restProc < 2) break;
        }
    }

    // Gets the total size of the sub-lattice(without faces)
    m_subLatticeSize = 1;
    for (int i = 0; i < 4; i++) {
        m_subLatticeSize *= m_N[i];
    }
    m_lattice = new Links[m_subLatticeSize];

    // If has a size of 2, we exit as that may produce poor results.
    for (int i = 0; i < 4; i++) {
        if (m_N[i] <= 2) {
            if (m_processRank == 0) {
                cout << "Error: lattice size of 2 or less are not allowed: "; // Due to instabilities, possibly?
                for (int j = 0; j < 4; j++) {
                    cout << m_N[j] << " ";
                }
                cout << " --> exiting."<< endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            exit(0);
        }
    }

    // Sets up number of processors per dimension
    for (int i = 0; i < 3; i++) {
        m_processorsPerDimension[i] = m_NSpatial / m_N[i];
    }
    m_processorsPerDimension[3] = m_NTemporal / m_N[3];

    // Initializes the neighbour lists
    m_neighbourLists->initialize(m_processRank, m_numprocs, m_processorsPerDimension);

    // Passes relevant information to the index handler(for the shifts).
    m_indexHandler->setN(m_N);
    m_indexHandler->setNTot(m_NSpatial, m_NTemporal);
    m_indexHandler->setNeighbourList(m_neighbourLists);

    // Passes the index handler and dimensionality to the action and correlator classes.
    m_S->initializeIndexHandler(m_indexHandler);
    m_S->setN(m_N);
    m_correlator->initializeIndexHandler(m_indexHandler);
    m_correlator->setN(m_N);
    m_correlator->setLatticeSize(m_subLatticeSize);

    // Volumes in sub dimension
    m_VSub[0] = m_N[0]; // X-dimension
    m_VSub[1] = m_N[0]*m_N[1]; // Y-dimension
    m_VSub[2] = m_N[0]*m_N[1]*m_N[2]; // Z-dimension
    m_VSub[3] = m_N[0]*m_N[1]*m_N[2]*m_N[3]; // T-dimension

    // Volumes in total lattice
    m_V[0] = m_NSpatial; // X volume
    m_V[1] = m_NSpatial*m_NSpatial; // XY volume
    m_V[2] = m_NSpatial*m_NSpatial*m_NSpatial; // XYZ Volume
    m_V[3] = m_NSpatial*m_NSpatial*m_NSpatial*m_NTemporal; // XYZT volume

    // TESTS ====================================================================================
    if (m_processRank == 0) {
        cout << m_numprocs << " number of threads" << endl;
        cout << "Processsors per dimension: ";
        for (int i = 0; i < 4; i++) {
            cout << m_processorsPerDimension[i] << " ";
        }
        cout << endl;
        cout << "Lattice dimensions: N = " << m_NSpatial << " N_T = " << m_NTemporal << endl;
        cout << "Lattice volumes: ";
        for (int i = 0; i < 4; i++) {
            cout << m_V[i] << " ";
        }
        cout << endl;
        cout << "Sub lattice dimensions: ";
        for (int i = 0; i < 4; i++) {
            cout << m_N[i] << " ";
        }
        cout << endl;
        cout << "Sub-lattice volumes: ";
        for (int i = 0; i < 4; i++) {
            cout << m_VSub[i] << " ";
        }
        cout << endl;
        cout << "Sub lattice size: " << m_subLatticeSize << endl;
        cout << "Link size: " << linkSize << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
//    m_neighbourLists->getNeighbours(m_processRank)->print();
//    cout << "Neighbourlist coordinates for P = " << m_processRank << " is: ";
//    for (int i = 0; i < 4; i++) {
//        cout << m_neighbourLists->getProcessorDimensionPosition(i) << " ";
//    }
//    cout << endl;
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank == 0) cout << "Exits at neighbour lists" << endl;
//    MPI_Finalize(); exit(0);
    // ==========================================================================================

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
        for (unsigned int x = 0; x < m_N[0]; x++) {
            for (unsigned int y = 0; y < m_N[1]; y++) {
                for (unsigned int z = 0; z < m_N[2]; z++) {
                    for (unsigned int t = 0; t < m_N[3]; t++) {
                        for (unsigned int mu = 0; mu < 4; mu++) {
//                            m_lattice[getIndex(x,y,z,t,m_N[1],m_N[2],m_N[3])].U[mu].identity();
                            m_lattice[m_indexHandler->getIndex(x,y,z,t)].U[mu].identity();
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
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice, x, y, z, t, mu);
                        for (int n = 0; n < m_NUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
//                            updateLink(getIndex(x, y, z, t, m_N[1], m_N[2], m_N[3]), mu);
                            updateLink(m_indexHandler->getIndex(x,y,z,t), mu);
                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu);
                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            {
//                                m_lattice[getIndex(x, y, z, t, m_N[1], m_N[2], m_N[3])].U[mu].copy(m_updatedMatrix);
                                m_lattice[m_indexHandler->getIndex(x,y,z,t)].U[mu].copy(m_updatedMatrix);
                                m_acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
}


void System::runMetropolis(bool storePreObservables, bool writeConfigsToFile)
{
    //// TESTS ==============================================================================
    MPI_Barrier(MPI_COMM_WORLD);
//    writeConfigurationToFile(0); // Writing initial config to file.
    loadFieldConfiguration("para8core061102.bin");
//    loadFieldConfiguration("para32core0.627734.bin");
    // "para8to32core0563387.bin"
//    writeConfigurationToFile(1);
//    loadFieldConfiguration("configs_profiling_run_beta6.000000_config1.bin");
    MPI_Barrier(MPI_COMM_WORLD);

    //// PRINTS LATTICES AT EDGES IN ORDER TO CHECK THAT WE HAVE CORRECT LATTICE ORDERING
//    if (m_processRank==0) { // CHECK THAT THE CORRECT LINKS ARE WHERE THE SHOULD BE
//        cout << "Rank " << m_processRank << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]-1,0,0,0)].U[0].print();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==0) { // CHECK THAT THE CORRECT LINKS ARE WHERE THE SHOULD BE
//        cout << "Rank " << m_processRank << " half of P=0 (should equal the P=1 end link for 32 procs." << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]/2-1,0,0,0)].U[0].print();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==1) {
//        cout << "Rank " << m_processRank << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]-1,0,0,0)].U[0].print();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==2) {
//        cout << "Rank " << m_processRank << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]-1,0,0,0)].U[0].print();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==1) {
//        cout << "Rank " << m_processRank << " half of P=1 (should equal the P=2 end link for 32 procs." << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]/2-1,0,0,0)].U[0].print();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==3) {
//        cout << "Rank " << m_processRank << " should equal last of P=1 for 16 procs." << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]-1,0,0,0)].U[0].print();
//    }

//    loadFieldConfiguration("unityScalar.bin");
//    loadFieldConfiguration("unity16cores.bin");
//    loadFieldConfiguration("unity32cores.bin");
//    loadFieldConfiguration("parallel32core16cube_plaquette0594052.bin"); // From parallel program version, 32 cores
//    loadFieldConfiguration("config32updated0612263.bin"); // From parallel program version, 32 cores
//    loadFieldConfiguration("parallel8core16cube16.bin"); // From parallel program version, 8 cores
//    loadFieldConfiguration("scalar16cubed16run1.bin"); // From scalar program version
    MPI_Barrier(MPI_COMM_WORLD);
    double corr = m_correlator->calculate(m_lattice);
    MPI_Allreduce(&corr, &corr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    corr /= double(m_numprocs);
    if (m_processRank == 0) cout << "Plaquette value: " << corr << endl << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    writeConfigurationToFile(1);
    loadFieldConfiguration("configs_profiling_run_beta6.000000_config1.bin");
    corr = m_correlator->calculate(m_lattice);
    MPI_Allreduce(&corr, &corr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    corr /= double(m_numprocs);
    if (m_processRank == 0) cout << "Plaquette value: " << corr << endl << endl;
    MPI_Finalize(); exit(1);
    //// ===================================================================================

    // Variables for checking performance of the update.
    clock_t preUpdate, postUpdate;
    double updateStorer = 0;

    // Calculating correlator before any updates have began.
    m_GammaPreThermalization[0] = m_correlator->calculate(m_lattice);

    // Summing and sharing correlator to all processors before any updates has begun
    MPI_Allreduce(&m_GammaPreThermalization[0], &m_GammaPreThermalization[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Dividing by the number of processors in order to get the correlator.
    m_GammaPreThermalization[0] /= double(m_numprocs);
    if (m_processRank == 0) cout << "Pre-thermialization correlator: " << m_GammaPreThermalization[0] << endl;

    // Running thermalization
    for (int i = 0; i < m_NTherm; i++)
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
        }
    }

    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
//    m_acceptanceCounter = double(m_acceptanceCounter)/double(m_numprocs);

    // Printing post-thermalization correlator and acceptance rate
    if (m_processRank == 0) {
        cout << "Termalization complete. Acceptance rate: " << double(m_acceptanceCounter)/double(4*m_latticeSize*m_NUpdates*m_NTherm) << endl;
        cout << "Post-thermialization correlator: " << m_GammaPreThermalization[m_NTherm] << endl;
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
        if (m_processRank == 0) cout << "Plaquette value: " << m_Gamma[alpha] << endl;
        if (writeConfigsToFile) writeConfigurationToFile(alpha);
        // TEST ============================================================
        MPI_Barrier(MPI_COMM_WORLD);
        loadFieldConfiguration("configs_profiling_run_beta6.000000_config0.bin");
        MPI_Barrier(MPI_COMM_WORLD);
        double corr = m_correlator->calculate(m_lattice);
        MPI_Allreduce(&corr, &corr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        corr /= double(m_numprocs);
        if (m_processRank == 0) cout << "Plaquette value: " << corr << endl << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize(); exit(1);
        // =================================================================
    }
    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    m_acceptanceCounter = double(m_acceptanceCounter)/double(m_numprocs);
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
        file << "beta" << m_beta << endl;
        file << "acceptanceCounter " << getAcceptanceRate() << endl;
        file << "NCor " << m_NCor << endl;
        file << "NCf " << m_NCf << endl;
        file << "NTherm " << m_NTherm << endl;
        file << "AverageGamma " << m_averagedGamma << endl;
        file << "VarianceGamma " << m_varianceGamma << endl;
        file << "stdGamma " << m_stdGamma << endl;
        if (preThermalizationGamma) {
            for (int i = 0; i < m_NTherm+1; i++) {
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
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_NUpdates*m_latticeSize*4); // Times 4 from the Lorentz indices
}

void System::writeConfigurationToFile(int configNumber)
{
    /*
     * C-method for writing out configuration to file.
     * Arguments:
     *  configNumber   : (int) configuration number
     */
    if (m_processRank == 0) cout << "Writing configuration number " << configNumber << " to file." << endl;

    MPI_File file;
    std::string filename = m_outputFolder + m_filename + "beta" + std::to_string(m_beta) + "_config" + std::to_string(configNumber) + ".bin";
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;


//    for (int t = 0; t < m_N[3]; t++) {
//        nt = (m_neighbourLists->getProcessorDimensionPosition(3) * m_N[3] + t);
//        for (int z = 0; z < m_N[2]; z++) {
//            nz = m_V[0] * (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z) + nt;
//            for (int y = 0; y < m_N[1]; y++) {
//                ny = m_V[1] * (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y) + nz;
//                for (int x = 0; x < m_N[0]; x++) {
//                    nx = m_V[2] * (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x) + ny;
//                    MPI_File_write_at(file, nx*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
//                }
//            }
//        }
//    }
    for (unsigned int t = 0; t < m_N[3]; t++) {
        nt = (m_neighbourLists->getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned int z = 0; z < m_N[2]; z++) {
//            nz = m_V[0] * (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z) + nt;
            nz = (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned int y = 0; y < m_N[1]; y++) {
//                ny = m_V[1] * (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y) + nz;
                ny = (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned int x = 0; x < m_N[0]; x++) {
//                    nx = m_V[2] * (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x) + ny;
//                    MPI_File_read_at(file, nx*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    nx = (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x);
                    MPI_File_write_at(file, m_indexHandler->getGlobalIndex(nx,ny,nz,nt)*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    MPI_File_close(&file);
    if (m_processRank == 0) cout << "Configuration number " << configNumber << " written to " << filename << endl;
}

void System::loadFieldConfiguration(std::string filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */

    bool writePythonFile = false;

    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, (m_outputFolder + filename).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;

    std::ofstream indexFile;
    if (writePythonFile) {
        indexFile.open("../python_scripts/parallel_" + std::to_string(m_numprocs) + "core_index_file_rank" + std::to_string(m_processRank) + ".dat");
    }

    int counter1 = 0;
    unsigned int counter2 = 0;
    unsigned long int counter3 = 0;
    for (unsigned int t = 0; t < m_N[3]; t++) {
        nt = (m_neighbourLists->getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned int z = 0; z < m_N[2]; z++) {
//            nz = m_V[0] * (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z) + nt;
            nz = (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned int y = 0; y < m_N[1]; y++) {
//                ny = m_V[1] * (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y) + nz;
                ny = (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned int x = 0; x < m_N[0]; x++) {
//                    nx = m_V[2] * (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x) + ny;
//                    MPI_File_read_at(file, nx*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    nx = (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x);

                    counter3 = m_indexHandler->getGlobalIndex(nx,ny,nz,nt);
                    counter2 = m_indexHandler->getGlobalIndex(nx,ny,nz,nt);
                    counter1 = m_indexHandler->getGlobalIndex(nx,ny,nz,nt);
                    if ((unsigned int) counter1 != counter2 || (unsigned long int) counter1 != counter3 || (unsigned long int) counter2 != counter3) {
                        cout << "INTEGER OVERFLOW" << endl;
                        cout << "Int: " << counter1 << " unsigned int: " << counter2 << "unsigned long int: " << counter3 << endl;
                    }

                    MPI_File_read_at(file, m_indexHandler->getGlobalIndex(nx,ny,nz,nt)*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    if (writePythonFile) {
                        indexFile << nx << " " << ny << " " << nz << " " << nt << " " << m_indexHandler->getGlobalIndex(nx,ny,nz,nt) << " " << m_processRank << endl;
                    }
                }
            }
        }
    }
//    if (m_processRank==0) { // Sanity check
//        cout << "First link: " << endl;
//        m_lattice[m_indexHandler->getIndex(0,0,0,0)].U[0].print();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==m_numprocs-1) { // Sanity check
//        cout << "Last link: " << endl;
//        m_lattice[m_indexHandler->getIndex(m_N[0]-1,m_N[1]-1,m_N[2]-1,m_N[3]-1)].U[3].print();
//    }

    if (writePythonFile) {
        indexFile.close();
        cout << "../python_scripts/parallel_" + std::to_string(m_numprocs) + "core_index_file_rank" + std::to_string(m_processRank) + ".dat" << " written" << endl;
    }

    MPI_File_close(&file);
    if (m_processRank == 0) cout << "Configuration " << m_outputFolder + filename << " loaded." << endl;

//    // IMPLEMENT LOAD FROM FILE FOR LATTICE HERE
//    FILE *file; // C method
//    file = fopen((m_inputFolder +"_p" + std::to_string(m_processRank) + filename).c_str(), "rb");
//    for (int t = 0; t < m_N[3]; t++) {
//        for (int z = 0; z < m_N[2]; z++) {
//            for (int y = 0; y < m_N[1]; y++) {
//                for (int x = 0; x < m_N[0]; x++) {
//                    for (int mu = 0; mu < 4; mu++) {
//                        fread(&m_lattice[getIndex(x, y, z, t, m_N[1], m_N[2], m_N[3])].U[mu],sizeof(SU3),1,file);
//                    }
//                }
//            }
//        }
//    }
//    fclose(file);
//    cout << m_inputFolder + filename  + " loaded" << endl;
}
