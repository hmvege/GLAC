/*!
 * \class Communicator
 *
 * \brief Communicator holds methods for performing single link communication, which is need in the Action.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
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
    static int m_processorsPerDimension[4];

    // Fetching variables
    static bool muDir;
    static bool nuDir;
    static SU3 exchangeU;

    // Sub lattice dimensions
    static std::vector<unsigned int> m_N;

    // Private fetchSU3 functions
    inline static void MPIfetchSU3Positive(Lattice<SU3> *lattice, std::vector<int> n, int mu, const int SU3Dir);
    inline static void MPIfetchSU3Negative(Lattice<SU3> *lattice, std::vector<int> n, int mu, const int SU3Dir);
public:
    Communicator();
    ~Communicator();

    // Initializers
    static void init(int *numberOfArguments, char ***cmdLineArguments);
    static void initializeSubLattice();

    // Link getters
    static SU3 getPositiveLink(Lattice<SU3> *lattice, std::vector<int> n, int mu, int *muIndex, const int SU3Dir);
    static SU3 getNegativeLink(Lattice<SU3> *lattice, std::vector<int> n, const int mu, int *muIndex, const int SU3Dir);
    static SU3 getNeighboursNeighbourLink(Lattice<SU3> * lattice, std::vector<int> n , int mu, int *muIndex, const int nu, int *nuIndex, const int SU3Dir);
    static SU3 getNeighboursNeighbourNegativeLink(Lattice<SU3> * lattice, std::vector<int> &n, const int mu, int *muIndex, const int nu, int *nuIndex, const int SU3Dir);

    // Getters
    static int getProcessRank() { return m_processRank; }
    static int getNumProc() { return m_numprocs; }

    // Setters
    static void setN(const std::vector<unsigned int> &N);

    // MPI
    static void MPIExit(const std::string &message);
    static void MPIPrint(const std::string &message);
    static void setBarrier();
    static void setBarrierActive();
    static void gatherDoubleResults(double * data, const unsigned int N);
    static void freeMPIGroups();

    // MPI method for reducing to a single dimension
    static void reduceToTemporalDimension(std::vector<double> &obsResults, const std::vector<double> &obs);

    // Validity checkers
    static void checkProcessorValidity();
    static void checkSubLatticeDimensionsValidity();
    static void checkSubLatticeValidity();

    // Function used for bug testing
    static void checkLattice(Lattice<SU3> *lattice, const std::string &message);
};
}

#endif // COMMUNICATOR_H
