#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include "math/lattice.h"

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
    inline static void MPIfetchSU3Positive(Lattice<SU3> *lattice, std::vector<int> n, int mu, int SU3Dir);
    inline static void MPIfetchSU3Negative(Lattice<SU3> *lattice, std::vector<int> n, int mu, int SU3Dir);
public:
    Communicator();
    ~Communicator();

    // Initializers
    static void init(int *numberOfArguments, char ***cmdLineArguments);
    static void initializeSubLattice();

    // Link getters
    static SU3 getPositiveLink(Lattice<SU3> *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
    static SU3 getNegativeLink(Lattice<SU3> *lattice, std::vector<int> n, int mu, int *muIndex, int SU3Dir);
    static SU3 getNeighboursNeighbourLink(Lattice<SU3> * lattice, std::vector<int> n , int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);
    static SU3 getNeighboursNeighbourNegativeLink(Lattice<SU3> * lattice, std::vector<int> n, int mu, int *muIndex, int nu, int *nuIndex, int SU3Dir);

    // Getters
    static int getProcessRank() { return m_processRank; }
    static int getNumProc() { return m_numprocs; }

    // Setters
    static void setN(std::vector<unsigned int> N);

    // MPI
    static void MPIExit(std::string message);
    static void MPIPrint(std::string message);
    static void setBarrier();
    static void gatherDoubleResults(double * data, int N);

//    // MPI Groups
//    static MPI_Group ACTIVE_GROUP,WORLD_GROUP;
//    static MPI_Comm ACTIVE_COMM;

//    // Variable for storing if processor is active(will always seek maximum possible of cores divisible to 2^n
//    static bool m_active;

    // Validity checkers
    static void checkProcessorValidity();
    static void checkSubLatticeDimensionsValidity();
    static void checkSubLatticeValidity();
};
}

#endif // COMMUNICATOR_H
