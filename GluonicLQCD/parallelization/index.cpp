#include "index.h"

#include <iostream>

unsigned int Parallel::Index::m_N[] = {};
unsigned int Parallel::Index::m_NTot[] = {};

Parallel::Index::Index()
{
    /*
     * Index initialiser. After initialisation, must set the int*N and m_neighbourLists manually.
     */
}

Parallel::Index::~Index()
{

}

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

