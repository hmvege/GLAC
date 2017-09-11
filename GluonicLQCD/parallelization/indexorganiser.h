#ifndef INDEXORGANISER_H
#define INDEXORGANISER_H

#include "neighbours.h"
#include "links.h"
#include "matrices/su3.h"
#include <vector>

class indexOrganiser
{
private:
    int m_numprocs;
    int m_processRank;
    int m_N[4];
    Neighbours * m_neighbourLists = nullptr;
    SU3 exchangeU;
public:
    indexOrganiser(int numprocs, int processRank);
    ~indexOrganiser();

    // Index getter
    int getIndex(int i, int j, int k, int l);

    // Link getters
    SU3 getPositiveLink(Links *lattice, std::__1::vector<bool, _Allocator> indexes, int *mu);
    Links getNegativeLink(Links *lattice, std::__1::vector<bool, _Allocator> indexes, int *mu);
    Links getNeighboursNeighbourLink(Links * lattice, std::__1::vector<bool, _Allocator> indexes , int *mu, int *nu);

    // Setters
    void setN(int *N);
    void setNeighbourList(Neighbours *neighbourLists) { m_neighbourLists = neighbourLists; }
};

#endif // INDEXORGANISER_H
