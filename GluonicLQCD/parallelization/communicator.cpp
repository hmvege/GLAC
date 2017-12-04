#include "communicator.h"
#include "neighbours.h"
#include "index.h"
#include "config/parameters.h"
#include <mpi.h>

// Internal variables
bool Parallel::Communicator::muDir = 0;
bool Parallel::Communicator::nuDir = 0;
SU3 Parallel::Communicator::exchangeU; // Carefull! This might give a bug!
std::vector<unsigned int> Parallel::Communicator::m_N = {0,0,0,0};
// Variables used externally
int Parallel::Communicator::m_processRank = 0;
int Parallel::Communicator::m_numprocs = 0;

using Parallel::Neighbours;

Parallel::Communicator::Communicator()
{
}

Parallel::Communicator::~Communicator()
{
}

void Parallel::Communicator::MPIfetchSU3Positive(Lattice<SU3> *lattice, std::vector<int> n, int mu, int SU3Dir)
{
    /*
     * Performs an MPI call to retrieve matrix in the positive direction.
     * Arguments:
     *  lattice     : the entire lattice passed
     *  n           : position vector
     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
     *  SU3Dir      : SU3 matrix direction at link
     */
    MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0],n[1],n[2],n[3])],18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*mu],0, // Send
            &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*mu+1],0,                                               // Receive
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void Parallel::Communicator::MPIfetchSU3Negative(Lattice<SU3> *lattice, std::vector<int> n, int mu, int SU3Dir)
{
    /*
     * Performs an MPI call to retrieve matrix in the negative direction.
     * Arguments:
     *  lattice     : the entire lattice passed
     *  n           : position vector
     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
     *  SU3Dir      : SU3 matrix direction at link
     */
    MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0],n[1],n[2],n[3])],18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*mu+1],0,  // Send
            &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*mu],0,                                            // Receive
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

SU3 Parallel::Communicator::getPositiveLink(Lattice<SU3> *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
{
    /*
     * Function for retrieving link in positive direction.
     * Takes:
     *  lattice     : constisting of links
     *  n           : position vector in lattice
     *  mu          : lorentz index mu
     *  muIndex     : lorentz "vector"
     *  SU3Dir      : which of the four SU3 matrices which we need
     */
    if ((n[mu]+muIndex[mu]) % m_N[mu] == 0) {
        n[mu] = 0;
        MPIfetchSU3Positive(lattice,n,mu,SU3Dir);
        return exchangeU;
    }
    else {
        return lattice[SU3Dir][Index::getIndex(n[0]+muIndex[0], n[1]+muIndex[1], n[2]+muIndex[2], n[3]+muIndex[3])];
    }
}

SU3 Parallel::Communicator::getNegativeLink(Lattice<SU3> *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
{
    /*
     * Function for retrieving link in negative direction.
     * Takes:
     *  lattice     : constisting of links
     *  n           : position vector in lattice
     *  mu          : lorentz index mu
     *  muIndex     : lorentz "vector"
     *  SU3Dir      : which of the four SU3 matrices which we need
     */
    if ((n[mu] - muIndex[mu] + m_N[mu]) % m_N[mu] == (m_N[mu] - 1)) {
        n[mu] = m_N[mu] - 1;
        MPIfetchSU3Negative(lattice,n,mu,SU3Dir);
        return exchangeU;
    }
    else {
        return lattice[SU3Dir][Index::getIndex(n[0]-muIndex[0], n[1]-muIndex[1], n[2]-muIndex[2], n[3]-muIndex[3])];
    }
}

SU3 Parallel::Communicator::getNeighboursNeighbourLink(Lattice<SU3> * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
{
    /*
     * Gets the neighbours neighbour link.
     * mu: positive direction
     * nu: negative direction
     *
     * Takes:
     *  lattice     : constisting of links
     *  n           : position vector in lattice
     *  mu          : lorentz index mu
     *  muIndex     : lorentz "vector"
     *  nu          : lorentz index nu
     *  nuIndex     : lorentz "vector"
     *  SU3Dir      : which of the four SU3 matrices which we need
     */
    muDir = (n[mu] + muIndex[mu]) % m_N[mu] == 0;
    nuDir = (n[nu] - nuIndex[nu] + m_N[nu]) % m_N[nu] == (m_N[nu] - 1);
    if (muDir && (!nuDir)) { // (muDir & ~nuDir)
        // Positive mu direction
        n[mu] = 0;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])],18,MPI_DOUBLE, Neighbours::getNeighbours(m_processRank)->list[2*mu],0,   // Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*mu+1],0,                                                                                              // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0]+muIndex[0],n[1]+muIndex[1],n[2]+muIndex[2],n[3]+muIndex[3])],18,MPI_DOUBLE, Neighbours::getNeighbours(m_processRank)->list[2*nu+1],0, // Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) { // muDir & nuDir
        // True edge case
        n[mu] = 0;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0],n[1],n[2],n[3])],18,MPI_DOUBLE, Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[SU3Dir][Index::getIndex(n[0]+muIndex[0]-nuIndex[0], n[1]+muIndex[1]-nuIndex[1], n[2]+muIndex[2]-nuIndex[2], n[3]+muIndex[3]-nuIndex[3])];
    }
}

SU3 Parallel::Communicator::getNeighboursNeighbourNegativeLink(Lattice<SU3> * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
{
    /*
     * Gets the neighbours neighbour link.
     * mu: negative direction
     * nu: negative direction
     *
     * Takes:
     *  lattice     : constisting of links
     *  n           : position vector in lattice
     *  mu          : lorentz index mu
     *  muIndex     : lorentz "vector"
     *  nu          : lorentz index nu
     *  nuIndex     : lorentz "vector"
     *  SU3Dir      : which of the four SU3 matrices which we need
     */
    muDir = (n[mu] - muIndex[mu] + m_N[mu]) % m_N[mu] == (m_N[mu] - 1);
    nuDir = (n[nu] - nuIndex[nu] + m_N[nu]) % m_N[nu] == (m_N[nu] - 1);
    if (muDir && (!nuDir)) { // (muDir & ~nuDir)
        // Positive mu direction
        n[mu] = m_N[mu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])],18,MPI_DOUBLE, Neighbours::getNeighbours(m_processRank)->list[2*mu],0,   // Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*mu+1],0,                                                                                               // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0]-muIndex[0],n[1]-muIndex[1],n[2]-muIndex[2],n[3]-muIndex[3])],18,MPI_DOUBLE, Neighbours::getNeighbours(m_processRank)->list[2*nu+1],0, // Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) { // muDir & nuDir
        // True edge case
        n[mu] = m_N[mu] - 1;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0],n[1],n[2],n[3])],18,MPI_DOUBLE, Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[SU3Dir][Index::getIndex(n[0]-muIndex[0]-nuIndex[0], n[1]-muIndex[1]-nuIndex[1], n[2]-muIndex[2]-nuIndex[2], n[3]-muIndex[3]-nuIndex[3])];
    }
}

void Parallel::Communicator::checkSubLatticeValidity()
{
    /*
     * Tests to ensures that the sub lattice is correctly divided.
     */
    bool latticeSizeError = false;
    for (int i = 0; i < 3; i++) {
        if (Parameters::getNSpatial() % m_N[i] != 0) {
            if (m_processRank == 0) cout << "Error: spatial dimension(s) is not correct";
            latticeSizeError = true;
        }
    }
    if (Parameters::getNTemporal() % m_N[3] != 0) {
        if (m_processRank == 0) cout << "Error: temporal dimension(s) is not correct";
        latticeSizeError = true;
    }
    if (Parameters::getSubLatticeSize()*m_numprocs != Parameters::getLatticeSize()) {
        if (m_processRank == 0) cout << "Error: volume of sub dimensions does not match the total lattice volume";
        latticeSizeError = true;
    }
    if (latticeSizeError) {
        if (m_processRank == 0) {
            std::string errMsg = "";
            errMsg += ": dimensions:  ";
            for (int j = 0; j < 4; j++) errMsg += std::to_string(m_N[j]) + " ";
            errMsg += " --> exiting.";
            MPIExit(errMsg);
        }
    }
}

void Parallel::Communicator::init(int *numberOfArguments, char ***cmdLineArguments)
{
    // Initializing parallelization, HIDE THIS?
    MPI_Init (numberOfArguments, cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &m_numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &m_processRank);

    if ((*numberOfArguments) != 2) {
        Parallel::Communicator::MPIExit("Error: please provide a json file to parse.");
    }
    checkProcessorValidity();
}

void Parallel::Communicator::checkProcessorValidity()
{
    /*
     * Exits if number of processors are odd.
     */
    if (m_numprocs % 2 != 0) {
        MPIExit("Error: odd number of processors --> exiting.");
    }
}

void Parallel::Communicator::checkSubLatticeDimensionsValidity()
{
    /*
     * Checks if the sublattice has any dimensions less than 2, as
     * that may cause instabilities with self-communications ect.
     */
    for (int i = 0; i < 4; i++) {
        if (m_N[i] <= 2) {
            if (m_processRank == 0) {
                std::string errMsg = "";
                errMsg += "Error: lattice size of 2 or less are not allowed: "; // Due to instabilities, possibly?
                for (int j = 0; j < 4; j++) errMsg += std::to_string(m_N[j]) + " ";
                errMsg += " --> exiting.";
                MPIExit(errMsg);
            }

        }
    }
}

void Parallel::Communicator::initializeSubLattice()
{
    int restProc = m_numprocs;
    int processorsPerDimension[4];
    double subLatticeSize = 1;
    // Only finds the sub lattice size iteratively if no preset value has been defined.
    if (!Parameters::getSubLatticePreset()) {
        // Sets up sub lattice dimensionality without any splitting
        for (int i = 0; i < 3; i++) {
            m_N[i] = Parameters::getNSpatial();
        }
        m_N[3] = Parameters::getNTemporal();
        // Iteratively finds and sets the sub-lattice dimensions
        while (restProc >= 2) {
            for (int i = 0; i < 4; i++) { // Counts from x to t
                m_N[i] /= 2;
                restProc /= 2;
                if (restProc < 2) break;
            }
        }
        // Sets the sub lattice dimensions in case they have not been set in the loading of the configuration.
        Parallel::Index::setN(m_N);
        Parameters::setN(m_N);
        setN(m_N);
    }
    // Gets the total size of the sub-lattice(without faces)
    for (int i = 0; i < 4; i++) {
        subLatticeSize *= m_N[i];
    }
    Parameters::setSubLatticeSize(subLatticeSize);
    // Ensures correct sub lattice dimensions
    checkSubLatticeValidity();
    // If has a size of 2, we exit as that may produce poor results.
    checkSubLatticeDimensionsValidity();
    // Sets up number of processors per dimension
    for (int i = 0; i < 3; i++) {
        processorsPerDimension[i] = Parameters::getNSpatial() / m_N[i];
    }
    processorsPerDimension[3] = Parameters::getNTemporal() / m_N[3];
    // Initializes the neighbour lists
    Parallel::Neighbours::initialize(m_processRank, m_numprocs, processorsPerDimension);
    Parameters::setProcessorsPerDimension(processorsPerDimension);
}

void Parallel::Communicator::setBarrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

void Parallel::Communicator::MPIExit(std::string message)
{
    if (m_processRank == 0) printf("\n%s", message.c_str());
    setBarrier();
    MPI_Finalize();
    exit(m_processRank);
}

void Parallel::Communicator::gatherDoubleResults(double * data, int N)
{
    double *tempData = new double[N]; // Possibly bad?! TEST
    MPI_Allreduce(data,tempData,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (int i = 0; i < N; i++) data[i] = tempData[i];
    delete [] tempData;
}

void Parallel::Communicator::setN(std::vector<unsigned int> N)
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}
