#include "index.h"

std::vector<unsigned int> Parallel::Index::m_N;
std::vector<long long> Parallel::Index::m_NTot;

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

void Parallel::Index::setNTot(unsigned int NSpatial, unsigned int NTemporal)
{
    /*
     * Function for setting the dimensionality of the total lattice.
     * Takes:
     *  N       : an array of ints of length 4, where each element is the dimension size of either x,y,z or t
     */
    m_NTot.resize(4);
    for (int i = 0; i < 3; i++) {
        m_NTot[i] = (long long) NSpatial;
    }
    m_NTot[3] = (long long) NTemporal;
}

