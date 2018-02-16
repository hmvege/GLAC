#include "communicator.h"
#include "neighbours.h"
#include "index.h"
#include "parallelparameters.h"
#include "config/parameters.h"
#include <mpi.h>
#include <cmath>

// Internal variables
bool Parallel::Communicator::muDir = 0;
bool Parallel::Communicator::nuDir = 0;
SU3 Parallel::Communicator::exchangeU; // Carefull! This might give a bug!
std::vector<unsigned int> Parallel::Communicator::m_N = {0,0,0,0};
// Variables used externally
int Parallel::Communicator::m_processRank = 0;
int Parallel::Communicator::m_numprocs = 0;
int Parallel::Communicator::m_processorsPerDimension[4] = {0,0,0,0};

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
            Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
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
            Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
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
                Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0]+muIndex[0],n[1]+muIndex[1],n[2]+muIndex[2],n[3]+muIndex[3])],18,MPI_DOUBLE, Neighbours::getNeighbours(m_processRank)->list[2*nu+1],0, // Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
                Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) { // muDir & nuDir
        // True edge case
        n[mu] = 0;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0],n[1],n[2],n[3])],18,MPI_DOUBLE, Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
                Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
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
                Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0]-muIndex[0],n[1]-muIndex[1],n[2]-muIndex[2],n[3]-muIndex[3])],18,MPI_DOUBLE, Neighbours::getNeighbours(m_processRank)->list[2*nu+1],0, // Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
                Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) { // muDir & nuDir
        // True edge case
        n[mu] = m_N[mu] - 1;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[SU3Dir][Index::getIndex(n[0],n[1],n[2],n[3])],18,MPI_DOUBLE, Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
                &exchangeU,18,MPI_DOUBLE,Neighbours::getNeighbours((Neighbours::getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
                Parallel::ParallelParameters::ACTIVE_COMM,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[SU3Dir][Index::getIndex(n[0]-muIndex[0]-nuIndex[0], n[1]-muIndex[1]-nuIndex[1], n[2]-muIndex[2]-nuIndex[2], n[3]-muIndex[3]-nuIndex[3])];
    }
}

void Parallel::Communicator::reduceToTemporalDimension(double * obsResults, double * obs)
{
    /*
     * Reduces flow results in matrixresults to a the temporal dimension
     */
    // Sets up the temporary buffer
    double temp[Parameters::getNTemporal()];
    for (int i = 0; i < Parameters::getNTemporal(); i++) {
        temp[i] = 0;
    }

    for (int iFlow = 0; iFlow < Parameters::getNFlows() + 1; iFlow++) {
        // Places obs values into temporary buffer
        for (unsigned int it = 0; it < m_N[3]; it++) {
            temp[it + m_N[3]*Neighbours::getProcessorDimensionPosition(3)] = obs[iFlow * m_N[3] + it];
        }

        // Reduces reduces results from temporary buffer to target obsResults
        MPI_Allreduce(temp, &obsResults[iFlow * Parameters::getNTemporal()], Parameters::getNTemporal(), MPI_DOUBLE, MPI_SUM, ParallelParameters::ACTIVE_COMM);

        // Resets temporary buffer
        for (int it = 0; it < Parameters::getNTemporal(); it++) {
            temp[it] = 0;
        }
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
        std::string errMsg = "";
        errMsg += ": dimensions:  ";
        for (int j = 0; j < 4; j++) errMsg += std::to_string(m_N[j]) + " ";
        errMsg += " --> exiting.";
        MPIExit(errMsg);
    }
}

void Parallel::Communicator::init(int *numberOfArguments, char ***cmdLineArguments)
{
    // Initializing parallelization, HIDE THIS?
    MPI_Init (numberOfArguments, cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &m_numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &m_processRank);

    if ((*numberOfArguments) != 2) {
        MPIExit("Error: please provide a json file to parse.");
    }

    int maxProcRank = m_numprocs;
    int binaryCounter = 1;
    int counter = 0;
    int rest = 0;

    while ((maxProcRank - binaryCounter) != 0) {
        binaryCounter = pow(2,counter);
        counter++;

        rest = maxProcRank % binaryCounter;
        if (rest != 0) {
            maxProcRank -= rest;
        }
    }

    if (m_processRank >= maxProcRank) {
        Parallel::ParallelParameters::active = false;
    } else {
        Parallel::ParallelParameters::active = true;
    }

    m_numprocs = maxProcRank;

    // Creating world group
    MPI_Comm_group(MPI_COMM_WORLD,&Parallel::ParallelParameters::WORLD_GROUP);

    // Create group based on active processors
    int activeProcs[m_numprocs];
    for (int i = 0; i < m_numprocs; i++) activeProcs[i] = i;
    MPI_Group_incl(Parallel::ParallelParameters::WORLD_GROUP,m_numprocs,activeProcs,&Parallel::ParallelParameters::ACTIVE_GROUP);

    // Creates a new communications group for all the active processors
    MPI_Comm_create_group(MPI_COMM_WORLD,Parallel::ParallelParameters::ACTIVE_GROUP,0,&Parallel::ParallelParameters::ACTIVE_COMM);

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
            std::string errMsg = "";
            errMsg += "Error: lattice size of 2 or less are not allowed: "; // Due to instabilities, possibly?
            for (int j = 0; j < 4; j++) errMsg += std::to_string(m_N[j]) + " ";
            errMsg += " --> exiting.";
            MPIExit(errMsg);
        }
    }
}

void Parallel::Communicator::initializeSubLattice()
{
    int restProc = m_numprocs;
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
        m_processorsPerDimension[i] = Parameters::getNSpatial() / m_N[i];
    }
    m_processorsPerDimension[3] = Parameters::getNTemporal() / m_N[3];

    // Initializes the neighbour lists
    Parallel::Neighbours::initialize(m_processRank, m_numprocs, m_processorsPerDimension);
    Parameters::setProcessorsPerDimension(m_processorsPerDimension);
}

void Parallel::Communicator::setBarrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

void Parallel::Communicator::freeMPIGroups()
{
    MPI_Group_free(&Parallel::ParallelParameters::WORLD_GROUP);
    MPI_Group_free(&Parallel::ParallelParameters::ACTIVE_GROUP);
    MPI_Comm_free(&Parallel::ParallelParameters::ACTIVE_COMM);
}

void Parallel::Communicator::MPIExit(std::string message)
{
    if (m_processRank == 0) printf("\n%s", message.c_str());
    freeMPIGroups();
    setBarrier();
    MPI_Finalize();
    exit(0);
}

void Parallel::Communicator::MPIPrint(std::string message)
{
    setBarrier();
    if (m_processRank == 0) printf("\n%s", message.c_str());
    setBarrier();
}

void Parallel::Communicator::gatherDoubleResults(double * data, int N)
{
    double tempData[N]; // Possibly bad?! TEST
    MPI_Allreduce(data,tempData,N,MPI_DOUBLE,MPI_SUM,Parallel::ParallelParameters::ACTIVE_COMM);
    for (int i = 0; i < N; i++) data[i] = tempData[i];
}

void Parallel::Communicator::setN(std::vector<unsigned int> N)
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}
