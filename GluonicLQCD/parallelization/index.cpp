#include "index.h"

Parallel::Index::Index()
{
    /*
     * Index initialiser. After initialisation, must set the int*N and m_neighbourLists manually.
     */
}

Parallel::Index::~Index()
{

}

void Parallel::Index::setN(std::vector<unsigned int> N)
{
    /*
     * Function for setting the dimensionality of the sublattice.
     * Takes:
     *  N       : an array of ints of length 4, where each element is the dimension size of either x,y,z or t
     */
    m_N = N;
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

