#include "neighbours.h"

Neighbours::Neighbours()
{
}

Neighbours::~Neighbours()
{
    std::cout << "DELETING NEIGHBOR LISTS..." << std::endl;
    delete [] m_neighbourLists;
}

void Neighbours::initialize(int processRank, int numproc, int * processorsPerDim) {
    m_processRank = processRank;
    m_numproc = numproc;
    m_Nx = processorsPerDim[0];
    m_Ny = processorsPerDim[1];
    m_Nz = processorsPerDim[2];
    m_Nt = processorsPerDim[3];
    m_neighbourLists = new NeighbourList[m_numproc];
    generateNeighbourList();
    m_P[0] = processRank % m_Nx;
    m_P[1] = (processRank / m_Nx) % m_Ny;
    m_P[2] = (processRank / (m_Nx * m_Ny)) % m_Nz;
    m_P[3] = (processRank / (m_Nx * m_Ny * m_Nz)) % m_Nt;
}

void Neighbours::generateNeighbourList()
{
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */
    for (int Np = 0; Np < m_numproc; Np++) {
        m_neighbourLists[Np].rank = Np;
        m_neighbourLists[Np].list[0] = getXMinusOne(Np);
        m_neighbourLists[Np].list[1] = getXPlusOne(Np);
        m_neighbourLists[Np].list[2] = getYMinusOne(Np);
        m_neighbourLists[Np].list[3] = getYPlusOne(Np);
        m_neighbourLists[Np].list[4] = getZMinusOne(Np);
        m_neighbourLists[Np].list[5] = getZPlusOne(Np);
        m_neighbourLists[Np].list[6] = getTMinusOne(Np);
        m_neighbourLists[Np].list[7] = getTPlusOne(Np);
    }
}

NeighbourList* Neighbours::getNeighbours(int Np) {
    /*
     * Must return a list. Therefore, returning a reference. Failure to do so, will result in segfault error.
     * Arguments:
     *  Np  : process rank
     */
    return &m_neighbourLists[Np];
}
