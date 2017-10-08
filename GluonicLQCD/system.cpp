#include <random>   // For Mersenne-Twister19937
#include <chrono>
//#include <ctime>
#include <cmath>    // For exp()
#include <fstream>
#include <iostream>
#include <iomanip>
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
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

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
        MPI_Finalize();
        exit(1);
    }
    else if (m_numprocs % 16 != 0) {
        if (m_processRank == 0) cout << "Warning: running for an number of prosessors not divisible by 16 may slow down performance due to the geometry of the lattice." << endl;
    }
    int restProc = m_numprocs;

    // Only finds the sub lattice size iteratively if no preset value has been defined.
    if (!m_subLatticeSizePreset) {
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
                for (int j = 0; j < 4; j++) cout << m_N[j] << " ";
                cout << " --> exiting."<< endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }

    // Checking if the sub lattice size acutally divide the lattice correctly.
    bool latticeSizeError = false;
    for (int i = 0; i < 3; i++) {
        if (m_NSpatial % m_N[i] != 0) {
            latticeSizeError = true;
        }
    }
    if (m_subLatticeSize*m_numprocs != m_latticeSize) {
        latticeSizeError = true;
    }
    if (m_NTemporal % m_N[3] != 0) {
        latticeSizeError = true;
    }
    if (latticeSizeError) {
        if (m_processRank == 0) {
            cout << "Error: sub-lattice size invalid: ";
            for (int j = 0; j < 4; j++) cout << m_N[j] << " ";
            cout << " --> exiting."<< endl;
        }
        MPI_Finalize();
        exit(0);
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
}

void System::setSubLatticeDimensions(int *NSub)
{
    /*
     * Function for specifying sub-lattice dimensions.
     * Arguments:
     *  (int*) NSub     : takes 4 integers, one integer for each sub-lattice dimension.
     */
    for (int i = 0; i < 4; i++) {
        m_N[i] = NSub[i];
    }
    m_subLatticeSizePreset = true;
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
                            m_lattice[m_indexHandler->getIndex(x,y,z,t)].U[mu].identity();
                        }
                    }
                }
            }
        }
    }
    if (m_processRank == 0) {
        cout << "\nLattice setup complete\n" << endl;
    }
}

void System::printRunInfo(bool verbose) {
    if (m_processRank == 0) {
        for (int i = 0; i < 60; i++) cout << "=";
        cout << endl;
        cout << "Threads:                               " << m_numprocs << endl;
        if (verbose) cout << "Lattice size:                          " << m_latticeSize << endl;
        cout << "Lattice dimensions(spatial, temporal): " << m_NSpatial << " " << m_NTemporal << endl;
        cout << "N configurations:                      " << m_NCf << endl;
        cout << "N correlation updates:                 " << m_NCor << endl;
        cout << "N thermalization updates:              " << m_NTherm << endl;
        cout << "N link updates:                        " << m_NUpdates << endl;
        cout << "Beta:                                  " << m_beta << endl;
        if (verbose) {
            cout << "SU3Eps:                                " << m_SU3Generator->getEpsilon() << endl;
            cout << "Sub lattice Size:                      " << m_subLatticeSize << endl;
            cout << "Sub latticedimensions:                 ";
            for (int i = 0; i < 4; i++) {
                cout << m_N[i] << " ";
            }
            cout << endl;
            cout << "Processsors per dimension:             ";
            for (int i = 0; i < 4; i++) {
                cout << m_processorsPerDimension[i] << " ";
            }
            cout << endl;
        }
        for (int i = 0; i < 60; i++) cout << "=";
        cout << endl;
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
//    m_updatedMatrix = m_SU3Generator->generateRandom()*m_lattice[latticeIndex].U[mu]; // Shorter method of updating matrix
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
                            updateLink(m_indexHandler->getIndex(x,y,z,t), mu);
//                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu);
//                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            if (exp(-m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu)) > m_uniform_distribution(m_generator))
                            {
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


void System::runMetropolis(bool storeThermalizationObservables, bool writeConfigsToFile)
{
    m_storeThermalizationObservables = storeThermalizationObservables;
    //// TESTS ==============================================================================
//    MPI_Barrier(MPI_COMM_WORLD);
    // Common files
    // loadFieldConfiguration("unityScalar.bin");
//     loadFieldConfiguration("unity16cores.bin");
//     loadFieldConfiguration("unity32cores.bin");
//    MPI_Barrier(MPI_COMM_WORLD);
//    double corr = m_correlator->calculate(m_lattice);
//    MPI_Allreduce(&corr, &corr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    corr /= double(m_numprocs);
//    if (m_processRank == 0) cout << "Plaquette value: " << corr << endl << endl;
//    MPI_Finalize(); exit(1);
    //// ===================================================================================
    if (m_processRank == 0) {
        cout << "Store thermalization observables:      ";
        if (m_storeThermalizationObservables) {
            cout << "TRUE" << endl;
        } else {
            cout << "FALSE" << endl;
        }
        cout << "Store configurations:                  ";
        if (writeConfigsToFile) {
            cout << "TRUE" << endl;
        } else {
            cout << "FALSE" << endl;
        }
        for (int i = 0; i < 60; i++) cout << "=";
        cout << endl;
    }

    // Variables for checking performance of the update.
    steady_clock::time_point preUpdate, postUpdate;
    duration<double> updateTime;
    double updateStorer = 0;

//    clock_t preUpdate, postUpdate;
    if (m_storeThermalizationObservables) {
        // Calculating correlator before any updates have began.
        m_GammaPreThermalization[0] = m_correlator->calculate(m_lattice);

        // Summing and sharing correlator to all processors before any updates has begun
        MPI_Allreduce(&m_GammaPreThermalization[0], &m_GammaPreThermalization[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Dividing by the number of processors in order to get the correlator.
        m_GammaPreThermalization[0] /= double(m_numprocs);
        if (m_processRank == 0) {
            printf("\ni    Plaquette   ");
            printf("\n%-4d %-12.8f",0,m_GammaPreThermalization[0]);
        }
    }

    // Running thermalization
    for (int i = 1; i < m_NTherm+1; i++)
    {
        // Pre timer
        preUpdate = steady_clock::now();

        update();

        // Post timer
        postUpdate = steady_clock::now();
        updateTime = duration_cast<duration<double>>(postUpdate-preUpdate);
        updateStorer += updateTime.count();
        if (i % 20 == 0) { // Avg. time per update every 10th update
            if (m_processRank == 0) {
                printf("\nAvgerage update time(every 10th): %f sec.", updateStorer/double(i));
            }
        }

        // Print correlator every somehting or store them all(useful when doing the thermalization).
        if (m_storeThermalizationObservables) {
            // Calculating the correlator
            m_GammaPreThermalization[i] = m_correlator->calculate(m_lattice);

            // Summing and sharing results across the processors
            MPI_Allreduce(&m_GammaPreThermalization[i], &m_GammaPreThermalization[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Turn off!

            // Averaging the results
            m_GammaPreThermalization[i] /= double(m_numprocs);
            if (m_processRank == 0) {
                printf("\n%-4d %-12.8f",i,m_GammaPreThermalization[i]);
            }
        }
    }

    // Taking the average of the acceptance rate across the processors.
    if (m_NTherm != 0) MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    // Printing post-thermalization correlator and acceptance rate
    if (m_processRank == 0) {
        if (m_NTherm != 0) printf("\nTermalization complete. Acceptance rate: %f",double(m_acceptanceCounter)/double(4*m_latticeSize*m_NUpdates*m_NTherm));
        printf("\ni    Plaquette    Avg.Update-time   Accept/reject");
    }

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;

    // For measuring main program time CPU TIME
    clock_t preUpdateMain;
    preUpdateMain = clock();
    double avgUpdateMain = 0;

    // Main part of algorithm
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            preUpdateMain = clock();
            update();
        }
        avgUpdateMain += double(clock()-preUpdateMain)/double(CLOCKS_PER_SEC);
        // Averaging the gamma values
        m_Gamma[alpha] = m_correlator->calculate(m_lattice);
        MPI_Allreduce(&m_Gamma[alpha], &m_Gamma[alpha], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        m_Gamma[alpha] /= double(m_numprocs);

        if (m_processRank == 0) {
            // Printing plaquette value
            printf("\n%-4d %-12.8f   %-15.8f",alpha,m_Gamma[alpha],avgUpdateMain/double((alpha+1)*m_NCor));
            // Adding the acceptance ratio
            if (alpha % 10 == 0) {
                printf(" %-13.8f", double(m_acceptanceCounter)/double(4*m_subLatticeSize*(alpha+1)*m_NUpdates*m_NCor));
            }
        }

        // Writing field config to file
        if (writeConfigsToFile) writeConfigurationToFile(alpha);
    }
    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    if (m_processRank == 0) {
        printf("\n");
        for (int i = 0; i < 60; i++) cout << "=";
        printf("\nSystem completed.");
        printf("\nAverage update time: %.6f sec.", avgUpdateMain/double(m_NCf*m_NCor));
        printf("\nTotal update time for %d updates: %.6f sec.\n", m_NCf*m_NCor, avgUpdateMain);
    }
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
    if (m_processRank == 0) {
        for (int i = 0; i < 60; i++) cout << "=";
        cout << endl;
        cout << "Average plaqutte:      " << m_averagedGamma << endl;
        cout << "Standard deviation:    " << m_stdGamma << endl;
        cout << "Variance:              " << m_varianceGamma << endl;
        for (int i = 0; i < 60; i++) cout << "=";
        cout << endl;
    }
}

void System::writeDataToFile(std::string filename)
{
    /*
     * For writing the raw Gamma data to file.
     * Arguments:
     *  filename                : to write stats to
     */
    if (m_processRank == 0) {
        std::ofstream file;
        std::string fname = m_pwd + m_outputFolder + filename + ".dat";
        file.open(fname);
        file << "beta " << m_beta << endl;
        file << "acceptanceCounter " << getAcceptanceRate() << endl;
        file << "NCor " << m_NCor << endl;
        file << "NCf " << m_NCf << endl;
        file << "NTherm " << m_NTherm << endl;
        file << std::setprecision(15) << "AverageGamma " << m_averagedGamma << endl;
        file << std::setprecision(15) << "VarianceGamma " << m_varianceGamma << endl;
        file << std::setprecision(15) << "stdGamma " << m_stdGamma << endl;
        if (m_storeThermalizationObservables) {
            for (int i = 0; i < m_NTherm+1; i++) {
                file << std::setprecision(15) << m_GammaPreThermalization[i] << endl;
            }
            file << endl;
        }
        for (int i = 0; i < m_NCf; i++) {
            file << std::setprecision(15) << m_Gamma[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
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

    MPI_File file;
    std::string filename = m_pwd + m_outputFolder + m_filename
                                            + "_beta" + std::to_string(m_beta)
                                            + "_spatial" + std::to_string(m_NSpatial)
                                            + "_temporal" + std::to_string(m_NTemporal)
                                            + "_threads" + std::to_string(m_numprocs)
                                            + "_config" + std::to_string(configNumber) + ".bin";

    // TEMP =====================================================================================================
//    if (m_processRank == 0) cout << "Writing filename: " << filename << endl;
    // ==========================================================================================================

    MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
//    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;

    for (unsigned int t = 0; t < m_N[3]; t++) {
        nt = (m_neighbourLists->getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned int z = 0; z < m_N[2]; z++) {
            nz = (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned int y = 0; y < m_N[1]; y++) {
                ny = (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned int x = 0; x < m_N[0]; x++) {
                    nx = (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x);
                    MPI_File_write_at(file, m_indexHandler->getGlobalIndex(nx,ny,nz,nt)*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                }
            }
        }
    }
    MPI_File_close(&file);
}

void System::loadFieldConfiguration(std::string filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */
    MPI_File file;
    MPI_File_open(MPI_COMM_SELF, (m_pwd + m_outputFolder + filename).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
//    MPI_File_open(MPI_COMM_WORLD, (m_pwd + m_outputFolder + filename).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;

    for (unsigned int t = 0; t < m_N[3]; t++) {
        nt = (m_neighbourLists->getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned int z = 0; z < m_N[2]; z++) {
            nz = (m_neighbourLists->getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned int y = 0; y < m_N[1]; y++) {
                ny = (m_neighbourLists->getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned int x = 0; x < m_N[0]; x++) {
                    nx = (m_neighbourLists->getProcessorDimensionPosition(0) * m_N[0] + x);
                    MPI_File_read_at(file, m_indexHandler->getGlobalIndex(nx,ny,nz,nt)*linkSize, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], linkDoubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    MPI_File_close(&file);
    if (m_processRank == 0) cout << "Configuration " << m_outputFolder + filename << " loaded." << endl;
}
