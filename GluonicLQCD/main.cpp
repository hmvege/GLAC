#include <iostream>
#include <ctime>
#include <chrono>
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
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

/*
 * TODO:
 * [ ] Enforce sub lattice cubes when possible(when allocating dimensions)
 * [ ] Add lattice cube sizes manually from cmd line
 * [ ] Add function for loading fields? Or make a seperate program? Should probably be done here.
 * [ ] (optional) Switch to CORRECT method syntax, foo --> m_foo
 * [ ] (optional) Check that the lattice is gauge invariant: M^-1 * U * M, see Gattinger intro on how to make gauge fields gauge invariant!
 * [ ] (optional) Add better test suites, one that prints FAIL if test fails!!
 */

int main(int numberOfArguments, char* cmdLineArguments[])
{
    // Initializing parallelization
    int numprocs, processRank;
    MPI_Init (&numberOfArguments, &cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

    // Constants by default initialization
    int N           = 16;            // Spatial lattice points.
    int N_T         = 32;            // Temporal lattice points.
    int NTherm      = 100;           // Thermalization.
    int NCor        = 10;            // Correlation updates.
    int NCf         = 20;           // Number of configurations to generate.
    int NUpdates    = 10;           // Number of link updates before moving on.
    double beta     = 6.0;          // Beta value(connected to the lattice spacing a).
    double SU3Eps   = 0.24;         // Epsilon spread for generating SU(3) matrices.
    std::string batchName               = "TEST_RUN_CONFIG";
    bool writeConfigsToFile             = true;
    bool storeThermalizationPlaquettes  = false;
    bool hotStart                       = false;
    double seed                         = std::time(nullptr) + double(processRank);                     // Random matrix seed. Defualt: current time.
    double metropolisSeed               = std::time(nullptr) + double(numprocs) + double(processRank);  // Metropolis seed. Defualt: current time.

    // Only needed if numberOfArguments > 13.
    int NSub[4];

    if (numberOfArguments > 1) {
        batchName = cmdLineArguments[1];
    }
    if (numberOfArguments > 2) { // Points for each lattice dimension.
        N           = atoi(cmdLineArguments[2]);
    }
    if (numberOfArguments > 3) { // Time dimension.
        N_T         = atoi(cmdLineArguments[3]);
    }
    if (numberOfArguments > 4) { // Number of times we will thermalize. In production runs this should be around 2000,
        NTherm      = atoi(cmdLineArguments[4]);
    }
    if (numberOfArguments > 5) { // Only keeping every 20th path. In production runs this will be 200.
        NCor        = atoi(cmdLineArguments[5]);
    }
    if (numberOfArguments > 6) { // Number of field configurations that will be generated. In production runs, this should be around 1000 for the smallest, 500 for the others.
        NCf         = atoi(cmdLineArguments[6]);
    }
    if (numberOfArguments > 7) { // Number of times a link will be updated. Because of the 10% acceptance ratio, 10 times means about 1 update.
        NUpdates    = atoi(cmdLineArguments[7]);
    }
    if (numberOfArguments > 8) { // Beta value, phenomelogically connected to lattice spacing a.
        beta        = atof(cmdLineArguments[8]);
    }
    if (numberOfArguments > 9) { // Epsilon used for generating SU(3) matrices
        SU3Eps      = atof(cmdLineArguments[9]);
    }
    if (numberOfArguments > 10) { // If we are to write the configurations to file. Default is TRUE.
        if (atoi(cmdLineArguments[10]) == 0) {
            writeConfigsToFile = false;
        }
    }
    if (numberOfArguments > 11) { // If we are to store the plaquettes from the thermalization. Default is TRUE.
        if (atoi(cmdLineArguments[11]) == 1) {
            storeThermalizationPlaquettes = true;
        }
    }
    if (numberOfArguments > 12) { // Hot start/cold start. Default is cold start.
        if (atoi(cmdLineArguments[12]) == 1) {
            hotStart = true;
        }
    }

    if (numberOfArguments > 13) {
        for (int i = 13; i < 17; i++) {
            NSub[i % 13] = atoi(cmdLineArguments[i]);
        }

    }

//    if (numberOfArguments > 13) { // RNG Seed.
//        seed = atof(cmdLineArguments[13]) + double(processRank);
//    }
//    if (numberOfArguments > 14) { // Metrolis seed.
//        metropolisSeed = atof(cmdLineArguments[14]) + double(numprocs) + double(processRank);
//    }

//    if (processRank == 0) {
//        testInverseMatrix(SU3Eps, seed, 1e3, false);
//        inversePerformanceTest(SU3Eps,seed,1e5);
//        runTestSuite();
//        MPI_Finalize();
//        exit(0);
//    }

//    int NSub[4] = {8,8,8,8};

    // Program timers
    steady_clock::time_point programStart, programEnd;
    duration<double> programTime;
    programStart = steady_clock::now();

    // Main program part
    SU3MatrixGenerator SU3Gen(SU3Eps, seed);
    Plaquette G;
    WilsonGaugeAction S(beta);
    System pureGauge(N, N_T, NCf, NCor, NTherm, NUpdates, beta, metropolisSeed, &G, &S, numprocs, processRank);
    if (numberOfArguments > 13) pureGauge.setSubLatticeDimensions(NSub);
    pureGauge.latticeSetup(&SU3Gen, hotStart);
    pureGauge.setConfigBatchName(batchName);
    pureGauge.printRunInfo(true); // Always print run info
    pureGauge.runMetropolis(storeThermalizationPlaquettes, writeConfigsToFile);
    pureGauge.runBasicStatistics();
    pureGauge.printAcceptanceRate();
    pureGauge.writeDataToFile(batchName);

    // Finalizing and printing time taken
    programEnd = steady_clock::now();
    programTime = duration_cast<duration<double>>(programEnd-programStart);
    if (processRank == 0) cout << "Program complete. Time used: " << programTime.count() << " sec" << endl;

    MPI_Finalize();
    return 0;
}
