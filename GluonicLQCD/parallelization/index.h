#ifndef INDEX_H
#define INDEX_H

#include "neighbours.h"
#include "math/links.h"
#include "math/matrices/su3.h"
#include <vector>

namespace Parallel {
    class Index
    {
    private:
        static bool muDir;
        static bool nuDir;
        static unsigned int m_N[4];
        static unsigned int m_NTot[4];
        static Neighbours * m_neighbourLists;
        static SU3 exchangeU;
        static int m_processRank;

        // Private fetchSU3 functions
        static void MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
        static void MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
    public:
        Index();
        ~Index();

        // Index getter
        static unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l);
        static unsigned int getGlobalIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l);

        // Link getters
        static SU3 getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
        static SU3 getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
        static SU3 getNeighboursNeighbourLink(Links * lattice, std::vector<int> n , int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);
        static SU3 getNeighboursNeighbourNegativeLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);

        // Setters
        static void setProcessRank(int processRank) { m_processRank = processRank; }
        static void setN(unsigned int *N);
        static void setNTot(int NSpatial, int NTemporal);
        static void setNeighbourList(Neighbours *neighbourLists) { m_neighbourLists = neighbourLists; }

        static void setBarrier();
    };
}
#endif // INDEX_H
