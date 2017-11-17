#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include <mpi.h>
#include "index.h"
#include "neighbours.h"
#include "math/latticemath.h"

#include "config/parameters.h"

namespace Parallel {
    class Communicator
    {
    private:
        // Parallel variables
        static int m_processRank;
        static int m_numprocs;

        // Fetching variables
        static bool muDir;
        static bool nuDir;
        static SU3 exchangeU;

        // Sub lattice dimensions
        static unsigned int m_N[4];

        // Private fetchSU3 functions
        static void MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
        static void MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
    public:
        Communicator();
        ~Communicator();

        // Public variables
        static Neighbours * m_neighbourLists;

        // Link getters
        static SU3 getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
        static SU3 getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
        static SU3 getNeighboursNeighbourLink(Links * lattice, std::vector<int> n , int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);
        static SU3 getNeighboursNeighbourNegativeLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);

        // Setters
        static void init(int numprocs, int processRank) { setNumproc(numprocs); setProcessRank(processRank); }
        static void setProcessRank(int processRank) { m_processRank = processRank; }
        static void setNumproc(int numprocs) { m_numprocs = numprocs; }
        static void setN(unsigned int *N);
        static void setNeighbourList(Neighbours *neighbourLists) { m_neighbourLists = neighbourLists; }

        // Getters
        static int getProcessRank() { return m_processRank; }
        static int getNumProc() { return m_numprocs; }

        // MPI
        static void setBarrier();
        static void gatherDoubleResults(double * data, int N);

        // Validity checkers
        static void checkProcessorValidity();
        static void checkSubLatticeDimensionsValidity();
        static void checkSubLatticeValidity();
    };
}

#endif // COMMUNICATOR_H
