#include "indexorganiser.h"
#include <mpi.h>
#include "links.h"
#include "matrices/su3.h"
#include <vector>

indexOrganiser::indexOrganiser(int numprocs, int processRank) : m_numprocs(numprocs), m_processRank(processRank)
{
    // Shorthand should be n()
}

indexOrganiser::~indexOrganiser()
{

}

SU3 indexOrganiser::getPositiveLink(Links *lattice, std::vector indexes, int *muIndex, int mu)
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
    // IDEA: TAKE A VECTOR OF INDEXES AND MAKE IT ACCESSIBLE BY MU IN ORDER TO FIGURE OUT IF WE ARE AT AN EDGE OR NOT
    if (indeces[])

    if ((i+m_N[0]) % m_N[0] == 0) {
        // positive x direction, 2*i + 1
        // m_neighbourLists->getNeighbours(m_processRank)->list[2*mu[0] + 1]
        return MPI_Sendrecv(lattice[getIndex(0,j,k,l)].U[mu],18,MPI_DOUBLE,);
    }
    else if ((j+m_N[1]) % m_N[1] == 0) {
        // positive y direction
    }
    else if ((k+m_N[2]) % m_N[2] == 0) {
        // positive z direction
    }
    else if ((l+m_N[3]) % m_N[3] == 0) {
        // positive t direction
    }
    else {
        return lattice[getIndex(i,j,k,l)];
    }
//    return getIndex((i+N[0]) % N[0], (j+N[1]) % N[1], (k+N[2]) % N[2], (l+N[3]) % N[3], N[1], N[2], N[3]);
}

Links indexOrganiser::getNegativeLink(Links *lattice, std::vector indexes, int *muIndex, int mu)
{
    if ((i+m_N[0]) % m_N[0] == 0) {
        // positive x direction, 2*i + 1
        // m_neighbourLists->getNeighbours(m_processRank)->list[2*mu[0]]
        return MPI_Sendrecv(lattice[getIndex(0,j,k,l)].U[mu],18,MPI_DOUBLE,);
        MPI_Sendrecv(   m_lattice[getIndex(1,y,z,t)].U,
                        18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[0],0,
                        exchangeU,
                        18,MPI_DOUBLE,m_processRank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else if ((j+m_N[1]) % m_N[1] == 0) {
        // positive y direction
    }
    else if ((k+m_N[2]) % m_N[2] == 0) {
        // positive z direction
    }
    else if ((l+m_N[3]) % m_N[3] == 0) {
        // positive t direction
    }
    else {
        return lattice[getIndex(i,j,k,l)];
    }
}

Links indexOrganiser::getNeighboursNeighbourLink(Links * lattice, std::vector indexes, int *mu, int *nu, int mu, int nu)
{
    // m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu[0]]))->list[2*nu[0]]
    if (mu[0]) {
        if (nu[1]) {

        }
        else if (nu[2]) {

        }
        else if (nu[3]) {

        }
        else {

        }
    }
    else if (mu[1]) {

    }
    else if (mu[2]) {

    }
    else if (mu[3]) {

    }
//    if ((i+N[0]) % N[0] == 0) {
//        if ((j+N[1]) % N[1] == (N[1] - 1)) {

//        }
//        else if ((k+N[2]) % N[2] == (N[2] - 1)) {

//        }
//        else if ((l+N[3]) % N[3] == (N[3] - 1)) {

//        }
//        else {

//        }
//    }
//    else if ((j+N[1]) % N[1] == 0) {
//        if ((i+N[0]) % N[0] == (N[0] - 1)) {

//        }
//        else if ((k+N[2]) % N[2] == (N[2] - 1)) {

//        }
//        else if ((l+N[3]) % N[3] == (N[3] - 1)) {

//        }
//        else {

//        }
//    }
//    else if ((k+N[2]) % N[2] == 0) {
//        if ((i+N[0]) % N[0] == (N[0] - 1)) {

//        }
//        else if ((k+N[2]) % N[2] == (N[2] - 1)) {

//        }
//        else if ((l+N[3]) % N[3] == (N[3] - 1)) {

//        }
//        else {

//        }
//    }
//    else if ((l+N[3]) % N[3] == 0) {

//    }
    else {
        return lattice[getIndex(i,j,k,l)];
    }
}

int indexOrganiser::getIndex(int i, int j, int k, int l)
{
    return (m_N[3]*(m_N[2]*(m_N[1]*i + j) + k) + l);
}
