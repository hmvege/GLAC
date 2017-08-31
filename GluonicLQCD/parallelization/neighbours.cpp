#include "neighbours.h"
#include "neighbourlist.h"

Neighbours::Neighbours(int processRank, int numproc, int xDimProc, int yDimProc, int zDimProc, int tDimProc)
{
    m_processRank = processRank;
    m_numproc = numproc;
    m_xDimProc = xDimProc;
    m_yDimProc = yDimProc;
    m_zDimProc = zDimProc;
    m_tDimProc = tDimProc;
    neighbourLists = new NeighbourList[m_numproc];
    generateNeighbourList();
}

Neighbours::~Neighbours()
{
    delete [] neighbourLists;
}

void Neighbours::generateNeighbourList()
{

}

