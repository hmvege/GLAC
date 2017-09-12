#include <iostream>
#include <ctime>
#include <mpi.h>
#include "system.h"
#include "actions/action.h"
#include "actions/wilsongaugeaction.h"
#include "correlators/plaquette.h"
#include "matrices/su3matrixgenerator.h"
#include "parallelization/indexorganiser.h"

#include "unittests.h"

using std::cout;
using std::endl;

/*
 * TODO:
 * [x] Add plaquette correlator
 * [x] Make actions more general!! Aka, create a Wilson action
 * [x] Change to updating random matrices by X=RST
 * [x] Change to such that time dimension is 2N
 * [x] Add determinant for SU3 matrices
 * [x] Create method for saving lattice configuration
 * [x] Create method for loading lattice configuration
 * [x] Fix bug in matrices
 * [x] Find bug so that N != N_T works.
 * [x] Add write each Plaquette/observable to file-function.
 * [x] Finish SU3 basic properties unit testing such that I dont have to compare by hand
 * [x] Rename metropolis.cpp --> system.cpp
 * [x] Add shifting parallelization
 * [ ] Update functions for reading and writing sublattices
 * [ ] Switch to CORRECT method syntax, foo --> m_foo
 * [ ] Check that the lattice is gauge invariant: M^-1 * U * M, see Gattinger intro on how to make gauge fields gauge invariant!
 * [ ] Add better test suites, one that prints FAIL if test fails!!
 */

int main(int numberOfArguments, char* cmdLineArguments[])
{
    // Initializing parallelization
    int numprocs, processRank;
    MPI_Init (&numberOfArguments, &cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

    // Constants
    int N           = 8;            // Points for each lattice dimension, 8 points in time dimension
    int N_T         = 8;            // Time dimension
    int NTherm      = 2;           // Number of times we are to thermalize, that is NTherm * NCor. Should increase to 30 at least.
    int NCor        = 10;           // Only keeping every 20th path
    int NCf         = 10;          // Number of configurations to retrieve
    double beta     = 6;            // Should be
    double SU3Eps   = 0.24;         // Epsilon used for generating SU(3) matrices
    double seed     = std::time(nullptr) + double(processRank);
    double metropolisSeed = std::time(nullptr) + 1.0 + double(processRank);
    bool storeThermalizationPlaquettes = true;
    bool hotStart   = false;

//    if (processRank == 0) {
//        testInverseMatrix(SU3Eps, seed, 1e3, false);
//        inversePerformanceTest(SU3Eps,seed,1e5);
//        runTestSuite();
//        MPI_Finalize();
//        exit(0);
//    }

    clock_t programStart, programEnd;
    programStart = clock();
    SU3MatrixGenerator SU3Gen(SU3Eps, seed);
    Plaquette G;
    WilsonGaugeAction S(beta);
    System pureGauge(N, N_T, NCf, NCor, NTherm, metropolisSeed, &G, &S, numprocs, processRank);
    pureGauge.latticeSetup(&SU3Gen, hotStart);
    pureGauge.runMetropolis(storeThermalizationPlaquettes);
    pureGauge.runBasicStatistics();
    pureGauge.printAcceptanceRate();
//    pureGauge.writeConfigurationToFile("configs_profiling_run");
//    pureGauge.writeDataToFile("../output/correlatorData.dat");

    programEnd = clock();
    if (processRank == 0) cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    MPI_Finalize();
    cout << "MPI FINALIZE COMPLETE" << endl;
    return 0;
}
