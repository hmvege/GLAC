#ifndef INDEXORGANISER_H
#define INDEXORGANISER_H

#include "neighbours.h"

class indexOrganiser
{
private:
    int m_numprocs;
    int m_processRank;
    int m_N[4];

    Neighbours * m_neighbourLists = nullptr;
public:
    indexOrganiser(int numprocs, int processRank);
    ~indexOrganiser();

    // Setters
    void setN(int *N);
    void setNeighbourList(Neighbours *neighbourLists) { m_neighbourLists = neighbourLists; }
};

#endif // INDEXORGANISER_H
