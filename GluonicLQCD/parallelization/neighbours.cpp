#include "neighbours.h"

int Parallel::Neighbours::m_processRank;
int Parallel::Neighbours::m_numproc;
int Parallel::Neighbours::m_Nx = 0;
int Parallel::Neighbours::m_Ny = 0;
int Parallel::Neighbours::m_Nz = 0;
int Parallel::Neighbours::m_Nt = 0; // Prosessors per dimension
int Parallel::Neighbours::m_P[4]; // Prosessor coordinate
std::vector<NeighbourList> Parallel::Neighbours::m_neighbourLists = {0,0,0,0};

Parallel::Neighbours::Neighbours()
{
}

Parallel::Neighbours::~Neighbours()
{
}

void Parallel::Neighbours::initialize(int processRank, int numproc, int * processorsPerDim) {
    m_processRank = processRank;
    m_numproc = numproc;
    m_Nx = processorsPerDim[0];
    m_Ny = processorsPerDim[1];
    m_Nz = processorsPerDim[2];
    m_Nt = processorsPerDim[3];
    m_neighbourLists.resize(m_numproc);
    generateNeighbourList();
    m_P[0] = processRank % m_Nx;
    m_P[1] = (processRank / m_Nx) % m_Ny;
    m_P[2] = (processRank / (m_Nx * m_Ny)) % m_Nz;
    m_P[3] = (processRank / (m_Nx * m_Ny * m_Nz)) % m_Nt;
}

void Parallel::Neighbours::generateNeighbourList()
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

std::vector<NeighbourList> Parallel::Neighbours::getNeighbours(int Np) {
    /*
     * Must return a list. Therefore, returning a reference. Failure to do so, will result in segfault error.
     * Arguments:
     *  Np  : process rank
     */
    return &m_neighbourLists[Np];
}
