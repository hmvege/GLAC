#ifndef INDEXORGANISER_H
#define INDEXORGANISER_H

#include "neighbours.h"
#include "links.h"
#include "matrices/su3.h"
#include <vector>

class IndexOrganiser
{
private:
    int m_numprocs;
    int m_processRank;
    int m_N[4];
    Neighbours * m_neighbourLists = nullptr;
    SU3 exchangeU;

    // Private fetchSU3 functions
    void MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
    void MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
public:
    IndexOrganiser(int numprocs, int processRank, Neighbours *neighbourLists, int *N);
    ~IndexOrganiser();

    // Index getter
    int getIndex(int i, int j, int k, int l);

    // Link getters
    SU3 getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
    SU3 getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
    SU3 getNeighboursNeighbourLink(Links * lattice, std::vector<int> n , int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);

    // Setters
    void setN(int *N);
    void setNeighbourList(Neighbours *neighbourLists) { m_neighbourLists = neighbourLists; }
};

#endif // INDEXORGANISER_H
