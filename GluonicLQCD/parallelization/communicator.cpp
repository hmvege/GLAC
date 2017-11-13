#include "communicator.h"

bool Parallel::Communicator::muDir = 0;
bool Parallel::Communicator::nuDir = 0;
unsigned int Parallel::Communicator::m_N[4];
Neighbours * Parallel::Communicator::m_neighbourLists = nullptr;
SU3 Parallel::Communicator::exchangeU;
int Parallel::Communicator::m_processRank = 0;
int Parallel::Communicator::m_numprocs = 0;


Parallel::Communicator::Communicator()
{
}

Parallel::Communicator::~Communicator()
{
}

void Parallel::Communicator::MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir)
{
    /*
     * Performs an MPI call to retrieve matrix in the positive direction.
     * Arguments:
     *  lattice     : the entire lattice passed
     *  n           : position vector
     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
     *  SU3Dir      : SU3 matrix direction at link
     */
    MPI_Sendrecv(&lattice[Index::getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0, // Send
            &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                               // Receive
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void Parallel::Communicator::MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir)
{
    /*
     * Performs an MPI call to retrieve matrix in the negative direction.
     * Arguments:
     *  lattice     : the entire lattice passed
     *  n           : position vector
     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
     *  SU3Dir      : SU3 matrix direction at link
     */
    MPI_Sendrecv(&lattice[Index::getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,  // Send
            &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,                                            // Receive
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

SU3 Parallel::Communicator::getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
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
        return lattice[Index::getIndex(n[0]+muIndex[0], n[1]+muIndex[1], n[2]+muIndex[2], n[3]+muIndex[3])].U[SU3Dir];
    }
}

SU3 Parallel::Communicator::getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
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
        return lattice[Index::getIndex(n[0]-muIndex[0], n[1]-muIndex[1], n[2]-muIndex[2], n[3]-muIndex[3])].U[SU3Dir];
    }
}

SU3 Parallel::Communicator::getNeighboursNeighbourLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
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
        MPI_Sendrecv(&lattice[Index::getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,   // Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                                                                              // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[Index::getIndex(n[0]+muIndex[0],n[1]+muIndex[1],n[2]+muIndex[2],n[3]+muIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*nu+1],0, // Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) { // muDir & nuDir
        // True edge case
        n[mu] = 0;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[Index::getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[Index::getIndex(n[0]+muIndex[0]-nuIndex[0], n[1]+muIndex[1]-nuIndex[1], n[2]+muIndex[2]-nuIndex[2], n[3]+muIndex[3]-nuIndex[3])].U[SU3Dir];
    }
}

SU3 Parallel::Communicator::getNeighboursNeighbourNegativeLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
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
        MPI_Sendrecv(&lattice[Index::getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,   // Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                                                                              // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[Index::getIndex(n[0]-muIndex[0],n[1]-muIndex[1],n[2]-muIndex[2],n[3]-muIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*nu+1],0, // Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) { // muDir & nuDir
        // True edge case
        n[mu] = m_N[mu] - 1;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[Index::getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[Index::getIndex(n[0]-muIndex[0]-nuIndex[0], n[1]-muIndex[1]-nuIndex[1], n[2]-muIndex[2]-nuIndex[2], n[3]-muIndex[3]-nuIndex[3])].U[SU3Dir];
    }
}

void Parallel::Communicator::setBarrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

//void Parallel::Communicator::reduceDoubleArray(double *var, int N)
//{
//    MPI_Allreduce(&var, &var, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//}
