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
        static Neighbours m_neighbourLists;

        // Initializers
        static void init(int numprocs, int processRank);
        static void initializeSubLattice();

        // Link getters
        static SU3 getPositiveLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
        static SU3 getNegativeLink(Links *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
        static SU3 getNeighboursNeighbourLink(Links * lattice, std::vector<int> n , int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);
        static SU3 getNeighboursNeighbourNegativeLink(Links * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);

        // Getters
        static int getProcessRank() { return m_processRank; }
        static int getNumProc() { return m_numprocs; }

        // Setters
        static void setN(unsigned int *N);

        // MPI
        static void MPIExit(std::__1::string message);
        static void setBarrier();
        static void gatherDoubleResults(double * data, int N);

        // Validity checkers
        static void checkProcessorValidity();
        static void checkSubLatticeDimensionsValidity();
        static void checkSubLatticeValidity();
    };
}

#endif // COMMUNICATOR_H
