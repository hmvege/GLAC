#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "neighbourlist.h"

class Neighbours
{
private:
    int m_processRank;
    int m_numproc;
    int m_xDimProc;
    int m_yDimProc;
    int m_zDimProc;
    int m_tDimProc;
    void generateNeighbourList();
    NeighbourList * neighbourLists;
public:
    Neighbours(int processRank, int numproc, int xDimProc, int yDimProc, int zDimProc, int tDimProc);
    ~Neighbours();
    // Neighbour list
    int getNeighbour();

    inline int getAllNeighbourListsIndex(int rank) {
        return
    }
};

#endif // NEIGHBOURS_H
