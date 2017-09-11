#include "indexorganiser.h"
#include "links.h"
#include "neighbours.h"
#include "matrices/su3.h"
#include <mpi.h>
#include <vector>

indexOrganiser::indexOrganiser(int numprocs, int processRank, Neighbours neighbourLists) : m_numprocs(numprocs), m_processRank(processRank)
{
    setNeighbourList(&neighbourLists);
    // Shorthand should be n()
}

indexOrganiser::~indexOrganiser()
{

}

SU3 indexOrganiser::getPositiveLink(Links *lattice, std::vector<int> n, int *muIndex, int mu)
{
    /*
     * TWO CASES:
     * - Neighbour link
     * - Neighbours neighbor link, wont matter the order of the mu's
     * REQUIREMENTS:
     * - Need to pass the lorentz index, as that determines the directionality.
     * - The positive/negative direction determined by wetter or not the index is 0 (negative) or N (positive) after the modulus has been performed
     * IF i + N % N == 0 positive direction
     * IF i + N & N == N-1
     * IDEA: make (for single shift) one function for positive direction and one for negative direction
     */
    // IDEA: TAKE A VECTOR OF n AND MAKE IT ACCESSIBLE BY MU IN ORDER TO FIGURE OUT IF WE ARE AT AN EDGE OR NOT

    if (n[mu] % m_N[mu] == 0) {
        exchangeU = lattice[getIndex(n[0],n[1],n[2],n[3])].U[mu];
        MPI_Sendrecv(&exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,   // Send
                &lattice[getIndex(n[0],n[1],n[2],n[3])].U[mu],18,MPI_DOUBLE,m_processRank,0,                    // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[getIndex(n[0],n[1],n[2],n[3])].U[mu];
    }
//    if (muIndex[i]) {
//    }
//    for (int i = 0; i < 4; i++) {
//    }


//    return getIndex((i+N[0]) % N[0], (j+N[1]) % N[1], (k+N[2]) % N[2], (l+N[3]) % N[3], N[1], N[2], N[3]);
}

SU3 indexOrganiser::getNegativeLink(Links *lattice, std::vector<int> n, int *muIndex, int mu)
{
    if ((n[mu] + m_N[mu]) % m_N[mu] == 0) {
        exchangeU = lattice[getIndex(n[0],n[1],n[2],n[3])].U[mu];
        MPI_Sendrecv(&exchangeU,18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0, // Send
                &lattice[getIndex(n[0],n[1],n[2],n[3])].U[mu],18,MPI_DOUBLE,m_processRank,0,                 // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[getIndex(n[0],n[1],n[2],n[3])].U[mu];
    }
//    if ((i+m_N[0]) % m_N[0] == 0) {
//        // positive x direction, 2*i + 1
//        // m_neighbourLists->getNeighbours(m_processRank)->list[2*mu[0]]
//        return MPI_Sendrecv(lattice[getIndex(0,j,k,l)].U[mu],18,MPI_DOUBLE,);
//        MPI_Sendrecv(   m_lattice[getIndex(1,y,z,t)].U,
//                        18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[0],0,
//                        exchangeU,
//                        18,MPI_DOUBLE,m_processRank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//    }
//    else if ((j+m_N[1]) % m_N[1] == 0) {
//        // positive y direction
//    }
//    else if ((k+m_N[2]) % m_N[2] == 0) {
//        // positive z direction
//    }
//    else if ((l+m_N[3]) % m_N[3] == 0) {
//        // positive t direction
//    }
//    else {
//        return lattice[getIndex(i,j,k,l)];
//    }
}

SU3 indexOrganiser::getNeighboursNeighbourLink(Links * lattice, std::vector<int> n, int *muIndex, int *nuIndex, int mu, int nu)
{
    // Position of targeted neighbour
    m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu];
}

int indexOrganiser::getIndex(int i, int j, int k, int l)
{
    return (m_N[3]*(m_N[2]*(m_N[1]*i + j) + k) + l);
}
