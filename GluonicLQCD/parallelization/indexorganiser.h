#ifndef INDEXORGANISER_H
#define INDEXORGANISER_H

#include "neighbours.h"
#include "links.h"
#include "matrices/su3.h"
#include <vector>

class IndexOrganiser
{
private:
    bool muDir;
    bool nuDir;
    unsigned int m_N[4];
    unsigned int m_NTot[4];
    Neighbours * m_neighbourLists = nullptr;
    SU3 exchangeU;

    // Private fetchSU3 functions
    void MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
    void MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
public:
    int m_processRank;
    IndexOrganiser(int processRank);
    ~IndexOrganiser();

    // Index getter
    unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l);
    unsigned int getGlobalIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l);

    // Link getters
    SU3 getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
    SU3 getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
    SU3 getNeighboursNeighbourLink(Links * lattice, std::vector<int> n , int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);
    SU3 getNeighboursNeighbourNegativeLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);

    // Setters
    void setN(unsigned int *N);
    void setNTot(int NSpatial, int NTemporal);
    void setNeighbourList(Neighbours *neighbourLists) { m_neighbourLists = neighbourLists; }
};

#endif // INDEXORGANISER_H
