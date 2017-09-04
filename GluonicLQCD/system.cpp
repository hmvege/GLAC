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

    // For parallelization
    m_neighbourLists = new Neighbours;

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
        m_dim[i] = m_N;
    }
    m_dim[3] = m_N_T;
//    // TEST==========================================================
//    if (m_processRank == 0) {
//        cout << "Processor: " << m_processRank << endl;
//        for (int i = 0; i < 4; i++) {
//            cout << m_dim[i] << endl;
//        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    // ==============================================================

    // Iteratively finds and sets the sub-lattice cube sizes
    while (restProc >= 2) {
//        for (int i = 3; i > 0; i--) { // Counts from t to x
        for (int i = 0; i < 4; i++) { // Conts from x to t
            m_dim[i] /= 2;
            restProc /= 2;
            if (restProc < 2) break;
        }
    }
    m_subLatticeSize = 1;
    m_trueSubLatticeSize = 1;
    for (int i = 0; i < 4; i++) {
        m_subLatticeSize *= m_dim[i]; // Gets the total size of the sub-lattice(without faces)
        m_extendedDim[i] = m_dim[i] + 2; // Adds a face
        m_trueSubLatticeSize *= m_extendedDim[i]; // Gets the total size of the sub-lattice(with faces)
    }
    m_lattice = new Links[m_trueSubLatticeSize];
    // Passes around updated indexes for the sublattices
    m_S->setN(m_extendedDim);
    m_correlator->setN(m_extendedDim);
    m_correlator->setLatticeSize(m_subLatticeSize);

    // Sets up number of processors per dimension
    for (int i = 0; i < 3; i++) {
        m_processorsPerDimension[i] = m_N / m_dim[i];
    }
    m_processorsPerDimension[3] = m_N_T / m_dim[3];
    m_neighbourLists->initialize(m_processRank, m_numprocs, m_processorsPerDimension);

    // PRINTS ===========================================================
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank==0) {
//        cout << "Process rank: " << m_processRank << endl;
//        cout << "m_subLatticeSize = " << m_subLatticeSize << endl;
//        cout << "m_trueSubLatticeSize = " << m_trueSubLatticeSize << endl;
//        for (int i = 0; i < 4; i++) {
//            cout << "Dim: " << i << " processors: " << m_processorsPerDimension[i] << endl;
//        }
//        cout << "Processor: " << m_processRank << endl;
//        for (int i = 0; i < 4; i++) {
//            cout << "direction " << i << " | m_dim     = " << m_dim[i] << " | m_trueDim = " << m_extendedDim[i] << endl;
//        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    m_neighbourLists->getNeighbours(m_processRank)->print();
//    MPI_Barrier(MPI_COMM_WORLD);
//    exit(0);
//    // ==================================================================
}

void System::latticeSetup(SU3MatrixGenerator *SU3Generator, bool hotStart)
{
    /*
     * Sets up the lattice and its matrices.
     */
    subLatticeSetup();

    // Also, set up blocks to use!!
    m_SU3Generator = SU3Generator;
    if (hotStart) {
        // All starts with a completely random matrix.
        for (int i = 0; i < m_latticeSize; i++)
//        for (int i = 0; i < m_trueSubLatticeSize; i++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                m_lattice[i].U[mu] = m_SU3Generator->generateRandom();
//            m_lattice[i].U[mu] = m_SU3Generator->generateRST();
            }
        }
    } else {
        for (int x = 1; x < m_extendedDim[0]-1; x++) {
            for (int y = 1; y < m_extendedDim[1]-1; y++) {
                for (int z = 1; z < m_extendedDim[2]-1; z++) {
                    for (int t = 1; t < m_extendedDim[3]-1; t++) {
                        for (int mu = 0; mu < 4; mu++) {
                            m_lattice[getIndex(x,y,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[mu] = m_SU3Generator->generateIdentity();
                        }
                    }
                }
            }
        }
//        for (int i = 0; i < m_trueSubLatticeSize; i++)
//        {
//            for (int mu = 0; mu < 4; mu++)
//            {
//                // All starts with a completely random matrix. Observable should be 1.
//                m_lattice[i].U[mu] = m_SU3Generator->generateIdentity();
//            }
//        }
    }
    shareFaces();
    shareCorners();
    if (m_processRank == 0) {
        cout << "Lattice setup complete" << endl;
    }
}

void System::shareFaces()
{
    /*
     * Function for sharing faces, edges and vertexes of the hypercubes of the different processors.
     *
     * Neighbour list values defined as:
     *  0: x-1 | 1: x+1
     *  2: y-1 | 3: y+1
     *  4: z-1 | 5: z+1
     *  6: t-1 | 7: t+1
     */

    // Share x=0 (x-1 direction)
    for (int y = 1; y < m_extendedDim[1]-1; y++) {
        for (int z = 1; z < m_extendedDim[2]-1; z++) {
            for (int t = 1; t < m_extendedDim[3]-1; t++) {
                MPI_Sendrecv(   &m_lattice[getIndex(1,y,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[0],0,
                                &m_lattice[getIndex(m_extendedDim[0]-1,y,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[1],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

    // Share x=Nx (x+1 direction)
    for (int y = 1; y < m_extendedDim[1]-1; y++) {
        for (int z = 1; z < m_extendedDim[2]-1; z++) {
            for (int t = 1; t < m_extendedDim[3]-1; t++) {
                MPI_Sendrecv(   &m_lattice[getIndex(m_extendedDim[0]-2,y,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[1],0,
                                &m_lattice[getIndex(0,y,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[0],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

    // Share y=0  (y-1 direction)
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int z = 1; z < m_extendedDim[2]-1; z++) {
            for (int t = 1; t < m_extendedDim[3]-1; t++) {
                MPI_Sendrecv(   &m_lattice[getIndex(x,1,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0], // Send
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2],0,
                                &m_lattice[getIndex(x,m_extendedDim[1]-1,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0], // Recieve
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[3],0,
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

    // Share y=Ny (y+1 direction)
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int z = 1; z < m_extendedDim[2]-1; z++) {
            for (int t = 1; t < m_extendedDim[3]-1; t++) {
                MPI_Sendrecv(   &m_lattice[getIndex(x,m_extendedDim[1]-2,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[3],0,
                                &m_lattice[getIndex(x,0,z,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

    // Share z=0  (z-1 direction)
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int y = 1; y < m_extendedDim[1]-1; y++) {
            for (int t = 1; t < m_extendedDim[3]-1; t++) {
                MPI_Sendrecv(   &m_lattice[getIndex(x,y,1,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0], // send
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[4],0,
                                &m_lattice[getIndex(x,y,m_extendedDim[2]-1,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0], // receive
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[5],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

    // Share z=Nz (z+1 direction)
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int y = 1; y < m_extendedDim[1]-1; y++) {
            for (int t = 1; t < m_extendedDim[3]-1; t++) {
                MPI_Sendrecv(   &m_lattice[getIndex(x,y,m_extendedDim[2]-2,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[5],0,
                                &m_lattice[getIndex(x,y,0,t,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[4],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

    // Share t=0  (t-1 direction)
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int y = 1; y < m_extendedDim[1]-1; y++) {
            for (int z = 1; z < m_extendedDim[2]-1; z++) {
                MPI_Sendrecv(   &m_lattice[getIndex(x,y,z,1,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[6],0,
                                &m_lattice[getIndex(x,y,z,m_extendedDim[3]-1,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[7],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);                
            }
        }
    }

    // Share t=Nt (t+1 direction)
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int y = 1; y < m_extendedDim[1]-1; y++) {
            for (int z = 1; z < m_extendedDim[2]-1; z++) {
                MPI_Sendrecv(   &m_lattice[getIndex(x,y,z,m_extendedDim[3]-2,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[7],0,
                                &m_lattice[getIndex(x,y,z,0,m_extendedDim[1],m_extendedDim[2],m_extendedDim[3])].U[0].mat[0].z[0],
                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[6],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }
}

void System::shareCorners()
{
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << m_processRank << "   " << m_neighbourLists->getNeighbours(m_processRank)->list[0] << endl;
//    exit(0);
//    if (m_processRank == 0) {
//        cout<<"rank 0"<<endl;
//        cout << m_subLatticeSize<<endl;
//        cout << (m_extendedDim[1]-2) * (m_extendedDim[2]-2) * (m_extendedDim[3]-2) << endl;
//        cout << (m_dim[1]*m_dim[2]*m_dim[3]) << endl;
//    }
//    if (m_processRank==0) {
//        int counter = 0;
//        int n1,n2,n3,n4;
//        int N1,N2,N3,N4;
//        for (int i = 0; i < 8; i++) { // Running over adjacent cube in each direction
//            cout <<"neighbours of neighbour in direction" << i << ":"<<endl;
//            for (int j = 0; j < 8; j++) { // Running over faces of cube to share(everyone exept the direction of the cube)
//                if (i/2 != j/2) {
//                    cout << "("<<i<<"/2) % 2=" << (i/2)%2<<endl;
////                    for (int n1 = 0; n1 < m_extendedDim[]; n1++) {
////                        for (int n2 = 0; n2 < m_extendedDim[]; n2++) {
////                            MPI_Sendrecv(   &m_lattice[].U[0].mat[0].z[0],
////                                            72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_neighbourLists->getNeighbours(m_processRank)->list[i])->list[j],0,
////                                            &m_lattice[].U[0].mat[0].z[0],
////                                            72,MPI_DOUBLE,m_processRank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
////                        }
////                    }
////                    cout << "i=" << i << " i/2=" << i / 2 << " | j=" << j << " j/2=" << j/2 << endl;
//                    cout << m_neighbourLists->getNeighbours(m_neighbourLists->getNeighbours(m_processRank)->list[i])->list[j] << endl;
//                    counter++;
////                    cout << i / 2 << " | " << j/2 << endl;
//                }
//            }
//        }
//        cout << counter << endl;
//    }
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
    for (int x = 1; x < m_extendedDim[0]-1; x++) {
        for (int y = 1; y < m_extendedDim[1]-1; y++) {
            for (int z = 1; z < m_extendedDim[2]-1; z++) {
                for (int t = 1; t < m_extendedDim[3]-1; t++) {
                    for (int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice, x, y, z, t, mu);
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(getIndex(x, y, z, t, m_extendedDim[1], m_extendedDim[2], m_extendedDim[3]), mu);
                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu);
                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[getIndex(x, y, z, t, m_extendedDim[1], m_extendedDim[2], m_extendedDim[3])].U[mu].copy(m_updatedMatrix);
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
    clock_t preUpdate, postUpdate;

    m_GammaPreThermalization[0] = m_correlator->calculate(m_lattice);
    MPI_Allreduce(&m_GammaPreThermalization[0], &m_GammaPreThermalization[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    m_GammaPreThermalization[0] /= double(m_numprocs);
    if (m_processRank == 0) {
        cout << "Pre-thermialization correlator:  " << m_GammaPreThermalization[0] << endl;
    }

    // Running thermalization
    for (int i = 0; i < m_NTherm*m_NCor; i++)
    {
        preUpdate = clock();
        update();
        postUpdate = clock();
//        if (m_processRank == 0) cout << "Mid update:  " << ((postUpdate - preUpdate)/((double)CLOCKS_PER_SEC)) << endl;
        shareFaces();
        postUpdate = clock();
//        if (m_processRank == 0) cout << "Post update: " << ((postUpdate - preUpdate)/((double)CLOCKS_PER_SEC)) << endl;
        // Print correlator every somehting or store them all(useful when doing the thermalization)
        if (storePreObservables) {
            m_GammaPreThermalization[i+1] = m_correlator->calculate(m_lattice);
            MPI_Allreduce(&m_GammaPreThermalization[i+1], &m_GammaPreThermalization[i+1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            m_GammaPreThermalization[i+1] /= double(m_numprocs);
            if (m_processRank == 0) cout << m_GammaPreThermalization[i+1] << endl;
        } else if ((i+1) % 10 == 0) {
            m_GammaPreThermalization[int((i+1)/10.0)+1] = m_correlator->calculate(m_lattice);
            MPI_Allreduce(&m_GammaPreThermalization[int((i+1)/10.0)+1], &m_GammaPreThermalization[int((i+1)/10.0)+1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            m_GammaPreThermalization[int((i+1)/10.0)+1] /= double(m_numprocs);
        }
    }

    if (m_processRank == 0) {
        cout << "Post-thermialization correlator: " << m_GammaPreThermalization[m_NTherm*m_NCor + 1] << endl;
        cout << "Termalization complete. Acceptance rate: " << m_acceptanceCounter/double(4*m_latticeSize*m_nUpdates*m_NTherm*m_NCor) << endl;
    }

    delete [] m_GammaPreThermalization; // De-allocating as this is not needed anymore

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;
    // Main part of algorithm
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
            // SHARE FACES HERE!?!
            shareFaces();
        }
        m_Gamma[alpha] = m_correlator->calculate(m_lattice);
        MPI_Allreduce(&m_Gamma[alpha], &m_Gamma[alpha], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        m_Gamma[alpha] /= double(m_numprocs);
    }
    if (m_processRank == 0) {
        cout << "System completed." << endl;
        cout << "Correlator = " << m_Gamma[m_NCf-1] << endl;
    }
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
    file = fopen((m_outputFolder + "_p" + std::to_string(m_processRank) + filename).c_str(), "wb");
    for (int t = 1; t < m_extendedDim[3]-1; t++) {
        for (int z = 1; z < m_extendedDim[2]-1; z++) {
            for (int y = 1; y < m_extendedDim[1]-1; y++) {
                for (int x = 1; x < m_extendedDim[0]-1; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fwrite(&m_lattice[getIndex(x, y, z, t, m_extendedDim[1], m_extendedDim[2], m_extendedDim[3])].U[mu],sizeof(SU3),1,file);
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
    file = fopen((m_inputFolder +"_p" + std::to_string(m_processRank) + filename).c_str(), "rb");
    for (int t = 1; t < m_extendedDim[3]-1; t++) {
        for (int z = 1; z < m_extendedDim[2]-1; z++) {
            for (int y = 1; y < m_extendedDim[1]-1; y++) {
                for (int x = 1; x < m_extendedDim[0]-1; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fread(&m_lattice[getIndex(x, y, z, t, m_extendedDim[1], m_extendedDim[2], m_extendedDim[3])].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }
    fclose(file);
    cout << m_inputFolder + filename  + " loaded" << endl;
}
