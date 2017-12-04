#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include "neighbours.h"
#include "math/links.h"

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
        static std::vector<unsigned int> m_N;

        // Private fetchSU3 functions
        inline static void MPIfetchSU3Positive(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
        inline static void MPIfetchSU3Negative(Links *lattice, std::vector<int> n, int mu, int SU3Dir);
    public:
        Communicator();
        ~Communicator();

        // Public variables
        static Neighbours m_NLists;

        // Initializers
        static void init(int *numberOfArguments, char ***cmdLineArguments);
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
        static void setN(std::vector<unsigned int> N);

        // MPI
        static void MPIExit(std::string message); // INLINE HERE!!!
        static void setBarrier();
        static void gatherDoubleResults(double * data, int N);

        // Validity checkers
        static void checkProcessorValidity();
        static void checkSubLatticeDimensionsValidity();
        static void checkSubLatticeValidity();
    };
}

#endif // COMMUNICATOR_H
