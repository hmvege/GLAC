#include "index.h"
#include <mpi.h>

//bool Parallel::Index::muDir = 0;
//bool Parallel::Index::nuDir = 0;
//Neighbours * Parallel::Index::m_neighbourLists = nullptr;
//SU3 Parallel::Index::exchangeU;
//int Parallel::Index::m_processRank = 0;

unsigned int Parallel::Index::m_N[] = {};
unsigned int Parallel::Index::m_NTot[] = {};

Parallel::Index::Index()
{
    /*
     * Index initialiser. After initialisation, must set the int*N and m_neighbourLists manually.
     * Takes:
     *  processRank : process ID as given by MPI initialisation
     */
}

Parallel::Index::~Index()
{

}

//void Parallel::Index::MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir)
//{
//    /*
//     * Performs an MPI call to retrieve matrix in the positive direction.
//     * Arguments:
//     *  lattice     : the entire lattice passed
//     *  n           : position vector
//     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
//     *  SU3Dir      : SU3 matrix direction at link
//     */
//    MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0, // Send
//            &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                               // Receive
//            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//}

//void Parallel::Index::MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir)
//{
//    /*
//     * Performs an MPI call to retrieve matrix in the negative direction.
//     * Arguments:
//     *  lattice     : the entire lattice passed
//     *  n           : position vector
//     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
//     *  SU3Dir      : SU3 matrix direction at link
//     */
//    MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,  // Send
//            &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,                                            // Receive
//            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//}

//SU3 Parallel::Index::getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
//{
//    /*
//     * Function for retrieving link in positive direction.
//     * Takes:
//     *  lattice     : constisting of links
//     *  n           : position vector in lattice
//     *  mu          : lorentz index mu
//     *  muIndex     : lorentz "vector"
//     *  SU3Dir      : which of the four SU3 matrices which we need
//     */

//    if ((n[mu]+muIndex[mu]) % m_N[mu] == 0) {
//        n[mu] = 0;
//        MPIfetchSU3Positive(lattice,n,mu,SU3Dir);
//        return exchangeU;
//    }
//    else {
//        return lattice[getIndex(n[0]+muIndex[0], n[1]+muIndex[1], n[2]+muIndex[2], n[3]+muIndex[3])].U[SU3Dir];
//    }
//}

//SU3 Parallel::Index::getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
//{
//    /*
//     * Function for retrieving link in negative direction.
//     * Takes:
//     *  lattice     : constisting of links
//     *  n           : position vector in lattice
//     *  mu          : lorentz index mu
//     *  muIndex     : lorentz "vector"
//     *  SU3Dir      : which of the four SU3 matrices which we need
//     */
//    if ((n[mu] - muIndex[mu] + m_N[mu]) % m_N[mu] == (m_N[mu] - 1)) {
//        n[mu] = m_N[mu] - 1;
//        MPIfetchSU3Negative(lattice,n,mu,SU3Dir);
//        return exchangeU;
//    }
//    else {
//        return lattice[getIndex(n[0]-muIndex[0], n[1]-muIndex[1], n[2]-muIndex[2], n[3]-muIndex[3])].U[SU3Dir];
//    }
//}

//SU3 Parallel::Index::getNeighboursNeighbourLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
//{
//    /*
//     * Gets the neighbours neighbour link.
//     * mu: positive direction
//     * nu: negative direction
//     *
//     * Takes:
//     *  lattice     : constisting of links
//     *  n           : position vector in lattice
//     *  mu          : lorentz index mu
//     *  muIndex     : lorentz "vector"
//     *  nu          : lorentz index nu
//     *  nuIndex     : lorentz "vector"
//     *  SU3Dir      : which of the four SU3 matrices which we need
//     */
//    muDir = (n[mu] + muIndex[mu]) % m_N[mu] == 0;
//    nuDir = (n[nu] - nuIndex[nu] + m_N[nu]) % m_N[nu] == (m_N[nu] - 1);
//    if (muDir && (!nuDir)) { // (muDir & ~nuDir)
//        // Positive mu direction
//        n[mu] = 0;
//        MPI_Sendrecv(&lattice[getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,   // Send
//                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                                                                              // Receive
//                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//        return exchangeU;
//    }
//    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
//        // Negative nu direction
//        n[nu] = m_N[nu] - 1;
//        MPI_Sendrecv(&lattice[getIndex(n[0]+muIndex[0],n[1]+muIndex[1],n[2]+muIndex[2],n[3]+muIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*nu+1],0, // Send
//                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
//                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//        return exchangeU;
//    }
//    else if (muDir && nuDir) { // muDir & nuDir
//        // True edge case
//        n[mu] = 0;
//        n[nu] = m_N[nu] - 1;
//        MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
//                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
//                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//        return exchangeU;
//    }
//    else {
//        return lattice[getIndex(n[0]+muIndex[0]-nuIndex[0], n[1]+muIndex[1]-nuIndex[1], n[2]+muIndex[2]-nuIndex[2], n[3]+muIndex[3]-nuIndex[3])].U[SU3Dir];
//    }
//}

//SU3 Parallel::Index::getNeighboursNeighbourNegativeLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
//{
//    /*
//     * Gets the neighbours neighbour link.
//     * mu: negative direction
//     * nu: negative direction
//     *
//     * Takes:
//     *  lattice     : constisting of links
//     *  n           : position vector in lattice
//     *  mu          : lorentz index mu
//     *  muIndex     : lorentz "vector"
//     *  nu          : lorentz index nu
//     *  nuIndex     : lorentz "vector"
//     *  SU3Dir      : which of the four SU3 matrices which we need
//     */
//    muDir = (n[mu] - muIndex[mu] + m_N[mu]) % m_N[mu] == (m_N[mu] - 1);
//    nuDir = (n[nu] - nuIndex[nu] + m_N[nu]) % m_N[nu] == (m_N[nu] - 1);
//    if (muDir && (!nuDir)) { // (muDir & ~nuDir)
//        // Positive mu direction
//        n[mu] = m_N[mu] - 1;
//        MPI_Sendrecv(&lattice[getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,   // Send
//                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                                                                              // Receive
//                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//        return exchangeU;
//    }
//    else if (nuDir && (!muDir)) { // (nuDir & ~muDir)
//        // Negative nu direction
//        n[nu] = m_N[nu] - 1;
//        MPI_Sendrecv(&lattice[getIndex(n[0]-muIndex[0],n[1]-muIndex[1],n[2]-muIndex[2],n[3]-muIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*nu+1],0, // Send
//                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*nu],0,                                                                                        // Receive
//                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//        return exchangeU;
//    }
//    else if (muDir && nuDir) { // muDir & nuDir
//        // True edge case
//        n[mu] = m_N[mu] - 1;
//        n[nu] = m_N[nu] - 1;
//        MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,// Send
//                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,                                             // Receive
//                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//        return exchangeU;
//    }
//    else {
//        return lattice[getIndex(n[0]-muIndex[0]-nuIndex[0], n[1]-muIndex[1]-nuIndex[1], n[2]-muIndex[2]-nuIndex[2], n[3]-muIndex[3]-nuIndex[3])].U[SU3Dir];
//    }
//}

unsigned int Parallel::Index::getIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    /*
     * Function for retrieving lattice position in contigious memory allocation.
     *  i   : x position
     *  j   : y position
     *  k   : z position
     *  l   : t position
     */
    return i + m_N[0]*(j + m_N[1]*(k + m_N[2]*l)); // column-major
}

unsigned int Parallel::Index::getGlobalIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    /*
     * Function for retrieving global lattice position.
     *  i   : x position
     *  j   : y position
     *  k   : z position
     *  l   : t position
     */
    return i + m_NTot[0]*(j + m_NTot[1]*(k + m_NTot[2]*l)); // column-major
}

void Parallel::Index::setN(unsigned int *N)
{
    /*
     * Function for setting the dimensionality of the sublattice.
     * Takes:
     *  N       : an array of ints of length 4, where each element is the dimension size of either x,y,z or t
     */
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}

void Parallel::Index::setNTot(int NSpatial, int NTemporal)
{
    /*
     * Function for setting the dimensionality of the total lattice.
     * Takes:
     *  N       : an array of ints of length 4, where each element is the dimension size of either x,y,z or t
     */
    for (int i = 0; i < 3; i++) {
        m_NTot[i] = (unsigned int) NSpatial;
    }
    m_NTot[3] = (unsigned int) NTemporal;
}
