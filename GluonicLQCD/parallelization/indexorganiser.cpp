#include "indexorganiser.h"
#include "links.h"
#include "neighbours.h"
#include "matrices/su3.h"
#include <mpi.h>
#include <vector>

// TEMP TESTS
#include <iostream>
using std::cout;
using std::endl;

IndexOrganiser::IndexOrganiser(int processRank) : m_processRank(processRank)
{
    /*
     * IndexOrganiser initialiser. After initialisation, must set the int*N and m_neighbourLists manually.
     * Takes:
     *  processRank : process ID as given by MPI initialisation
     */
}

IndexOrganiser::~IndexOrganiser()
{

}

void IndexOrganiser::MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir)
{
    /*
     * Performs an MPI call to retrieve matrix in the positive direction.
     * Arguments:
     *  lattice     : the entire lattice passed
     *  n           : position vector
     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
     *  SU3Dir      : SU3 matrix direction at link
     */
    MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0, // Send
            &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                               // Receive
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void IndexOrganiser::MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir)
{
    /*
     * Performs an MPI call to retrieve matrix in the negative direction.
     * Arguments:
     *  lattice     : the entire lattice passed
     *  n           : position vector
     *  mu          : lorentz index for shift direction(always negative in either x,y,z or t direction)
     *  SU3Dir      : SU3 matrix direction at link
     */
    MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,  // Send
            &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,                                            // Receive
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

SU3 IndexOrganiser::getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
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

    //    if (mu < 0 || mu > 3) {
//        cout << "ERROR" << endl;
//        exit(1);
//    }
//    for (int i = 0; i < 4; i++) {
//        if (n[i] < 0 || n[i] > m_N[i]) {
//            cout << "ERROR IN n" << endl;
//            exit(1);
//        }
//    }

    if ((n[mu]+muIndex[mu]) % m_N[mu] == 0) {
        n[mu] = 0;

        //// ARE WE SHARING CORRECT?! --> Yes we are, apperantely. But can we make a test to ensure this is true?

        // Process rank sending to N+1 and receiving from N-1.
//        MPI_Barrier(MPI_COMM_WORLD);
//        cout << m_processRank << " " << 2*mu+1 << " " << m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1] << " " << 2*mu << " " << m_neighbourLists->getNeighbours(m_processRank)->list[2*mu] << endl;

//        MPI_Barrier(MPI_COMM_WORLD);
//        if (m_processRank==0) {
//            for (int i = 0; i < 4; i++) cout << n[0];
//            cout << "Process 0, positive direction " << "mu = " << 2*SU3Dir+1 << ", before sharing with processor " << m_neighbourLists->getNeighbours(0)->list[2*mu+1] << endl;
//            exchangeU.print();
//        }
//        int p=0;
//        double compareValue = 0;
//        int errCounter = 0;
//        if (m_processRank==m_neighbourLists->getNeighbours(p)->list[2*mu+1]) {
//            cout << "Process " << m_neighbourLists->getNeighbours(p)->list[2*mu+1] << " we are fetching SU3 matrix from for process 0." << endl;
//            lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir].print();
//        }

//        MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir].mat[0],1,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,&compareValue,1,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        /*
     MPI_Sendrecv(   m_lattice[getIndex(m_N[0]-2,y,z,t,m_N[1],m_N[2],m_N[3])].U,
    -                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[1],0,
    -                                m_lattice[getIndex(0,y,z,t,m_N[1],m_N[2],m_N[3])].U,
    -                                72,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[0],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    */

        MPIfetchSU3Positive(lattice,n,mu,SU3Dir);

//        if (SU3Dir==0 && m_processRank==0 && n[1] == 0 && n[2] == 0 && n[3] == 0) {
//            exchangeU.print();
//        }
        if (m_processRank) {
            if (exchangeU.mat[0].re() == 0.00719493) {
                cout << "Exchanged matrix at rank: " << m_processRank << endl;
                exchangeU.print();}
        }

//        MPI_Barrier(MPI_COMM_WORLD);
//        if (compareValue != lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir].mat[0].re()) {
//            errCounter++;
//            cout << compareValue << endl;
//            lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir].print();
//            exit(1);
//        }

//        if (errCounter != 0) {
//            cout << "Alert! Errors(?) detected: " << errCounter << endl;
//        }
//        exchangeU.identity(); // NOT PASSING CORRECTLY?
//        exchangeU.zeros();

//        MPI_Barrier(MPI_COMM_WORLD);
//        if (m_processRank==p) {
//            cout << "Process " << p << ", positive direction " << "mu = " << 2*SU3Dir+1 << ", after sharing with " << m_neighbourLists->getNeighbours(p)->list[2*mu+1] << endl;
//            exchangeU.print();
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        MPI_Finalize();exit(1);

        return exchangeU;
    }
    else {
        return lattice[getIndex(n[0]+muIndex[0], n[1]+muIndex[1], n[2]+muIndex[2], n[3]+muIndex[3])].U[SU3Dir];
    }
}

SU3 IndexOrganiser::getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir)
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
        return lattice[getIndex(n[0]-muIndex[0], n[1]-muIndex[1], n[2]-muIndex[2], n[3]-muIndex[3])].U[SU3Dir];
    }
}

SU3 IndexOrganiser::getNeighboursNeighbourLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir)
{
    /*
     * For our program, mu is always in the positive direction and nu is always in the negative direction.
     */
    bool muDir = (n[mu] + muIndex[mu]) % m_N[mu] == 0;
    bool nuDir = (n[nu] - nuIndex[nu] + m_N[nu]) % m_N[nu] == (m_N[nu] - 1);
    if (muDir && (!nuDir)) {
        // Positive mu direction
        n[mu] = 0;
        MPI_Sendrecv(&lattice[getIndex(n[0]-nuIndex[0],n[1]-nuIndex[1],n[2]-nuIndex[2],n[3]-nuIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1],0,   // Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*mu],0,                                                                                              // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (nuDir && (!muDir)) {
        // Negative nu direction
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[getIndex(n[0]+muIndex[0],n[1]+muIndex[1],n[2]+muIndex[2],n[3]+muIndex[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours(m_processRank)->list[2*nu],0, // Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours(m_processRank)->list[2*nu+1],0,                                                                                        // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else if (muDir && nuDir) {
        // True edge case
        n[mu] = 0;
        n[nu] = m_N[nu] - 1;
        MPI_Sendrecv(&lattice[getIndex(n[0],n[1],n[2],n[3])].U[SU3Dir],18,MPI_DOUBLE, m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu+1]))->list[2*nu],0,// Send
                &exchangeU,18,MPI_DOUBLE,m_neighbourLists->getNeighbours((m_neighbourLists->getNeighbours(m_processRank)->list[2*mu]))->list[2*nu+1],0,                                             // Receive
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        return exchangeU;
    }
    else {
        return lattice[getIndex(n[0]+muIndex[0]-nuIndex[0], n[1]+muIndex[1]-nuIndex[1], n[2]+muIndex[2]-nuIndex[2], n[3]+muIndex[3]-nuIndex[3])].U[SU3Dir];
    }
}

unsigned int IndexOrganiser::getIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    /*
     * Function for retrieving lattice position in contigious memory allocation.
     *  i   : x position
     *  j   : y position
     *  k   : z position
     *  l   : t position
     */
//    return (m_N[3]*(m_N[2]*(m_N[1]*i + j) + k) + l); // row major
//    return l + m_N[3]*(k + m_N[2]*(j + m_N[1]*i)); // same as old, row major
    return i + m_N[0]*(j + m_N[1]*(k + m_N[2]*l)); // column-major
}

unsigned int IndexOrganiser::getGlobalIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    /*
     * Function for retrieving global lattice position.
     *  i   : x position
     *  j   : y position
     *  k   : z position
     *  l   : t position
     */
//    return (m_NTot[3]*(m_NTot[2]*(m_NTot[1]*i + j) + k) + l); // row-major
    return i + m_NTot[0]*(j + m_NTot[1]*(k + m_NTot[2]*l)); // column-major
}

void IndexOrganiser::setN(unsigned int *N)
{
    /*
     * Function for setting the dimensionality of the sublattice.
     * Takes:
     *  N       : an array of ints of length 4, where each element is the dimension size of either x,y,z or t
     */
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
    if (m_processRank == 0) {
        cout << "Setting up the sublattice(indexOrganizer.cpp): ";
        for (int i = 0; i < 4; i++) {
            cout << m_N[i] << " ";
        }
        cout << endl;
    }
}

void IndexOrganiser::setNTot(int NSpatial, int NTemporal)
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
