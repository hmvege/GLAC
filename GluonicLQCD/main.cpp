#include <iostream>
#include <ctime>
#include <chrono>
#include <mpi.h>
#include "system.h"
#include "actions/action.h"
#include "actions/wilsongaugeaction.h"
#include "observables/plaquette.h"
#include "math/matrices/su3matrixgenerator.h"
#include "parallelization/index.h"
#include "parameters/parameters.h"

#include "tests/unittests.h"
#include "tests/testsuite.h"

using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

/*
 * TODO:
 * [ ] Implement better map structure system. e.g. latmath.h ect
 * [ ] Enforce sub lattice cubes as standard
 * [ ] Add batch name to print-out
 * [ ] Add function for loading fields? Or make a seperate program? Should probably be done here.
 * [ ] (optional) Switch to CORRECT method syntax, foo --> m_foo
 * [ ] (optional) Check that the lattice is gauge invariant: M^-1 * U * M, see Gattinger intro on how to make gauge fields gauge invariant!
 */

void runUnitTests(SU3MatrixGenerator *SU3Gen, bool verbose);

int main(int numberOfArguments, char* cmdLineArguments[])
{
    // Initializing parallelization
    int numprocs, processRank;
    MPI_Init (&numberOfArguments, &cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

    // Constants by default initialization
    int N           = 16;            // Spatial lattice points.
    int N_T         = 16;            // Temporal lattice points.
    int NTherm      = 100;           // Thermalization.
    int NCor        = 10;            // Correlation updates.
    int NCf         = 20;           // Number of configurations to generate.
    int NUpdates    = 10;           // Number of link updates before moving on.
    int NFlows      = 0;            // Number of times to flow a variable.
    double beta     = 6.0;          // Beta value(connected to the lattice spacing a).
    double SU3Eps   = 0.24;         // Epsilon spread for generating SU(3) matrices.
    std::string batchName               = "TEST_RUN_CONFIG";
    std::string pwd                     = "/Users/hansmathiasmamenvege/Programming/FYSSP100/GluonAction"; // MAC
    bool writeConfigsToFile             = true;
    bool storeThermalizationPlaquettes  = false;
    bool hotStart                       = false;
    double seed                         = std::time(nullptr) + double(processRank);                     // Random matrix seed. Defualt: current time.
    double metropolisSeed               = std::time(nullptr) + double(numprocs) + double(processRank);  // Metropolis seed. Defualt: current time.
    bool runUnitTestsFlag               = false;
    bool runVerboseUnitTestsFlag        = false;

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
    if (numberOfArguments > 7) { // Number of flows per configuration
        NFlows      = atoi(cmdLineArguments[7]);
    }
    if (numberOfArguments > 8) { // Number of times a link will be updated. Because of the 10% acceptance ratio, 10 times means about 1 update.
        NUpdates    = atoi(cmdLineArguments[8]);
    }
    if (numberOfArguments > 9) { // Beta value, phenomelogically connected to lattice spacing a.
        beta        = atof(cmdLineArguments[9]);
    }
    if (numberOfArguments > 10) { // Epsilon used for generating SU(3) matrices
        SU3Eps      = atof(cmdLineArguments[10]);
    }
    if (numberOfArguments > 11) { // If we are to write the configurations to file. Default is TRUE.
        if (atoi(cmdLineArguments[11]) == 0) {
            writeConfigsToFile = false;
        }
    }
    if (numberOfArguments > 12) { // If we are to store the plaquettes from the thermalization. Default is TRUE.
        if (atoi(cmdLineArguments[12]) == 1) {
            storeThermalizationPlaquettes = true;
        }
    }
    if (numberOfArguments > 13) { // Hot start/cold start. Default is cold start.
        if (atoi(cmdLineArguments[13]) == 1) {
            hotStart = true;
        }
    }
    if (numberOfArguments > 14) { // Takes absolute file path
        pwd = cmdLineArguments[14];
    }
    if (numberOfArguments > 15) { // Takes specific sub lattice dimensions
        for (int i = 15; i < 19; i++) {
            NSub[i % 15] = atoi(cmdLineArguments[i]);
        }
    }
    if (numberOfArguments > 19) { // Flag for running unit tests
        if (atoi(cmdLineArguments[19]) == 1) {
            runUnitTestsFlag = true;
        }
    }
    if (numberOfArguments > 20) { // Flag for running unit tests
        if (atoi(cmdLineArguments[20]) == 1) {
            runVerboseUnitTestsFlag = true;
        }
    }

    // Program timers
    steady_clock::time_point programStart, programEnd;
    duration<double> programTime;
    programStart = steady_clock::now();

    // Parameter setting
    Parameters::setBatchName(batchName);
    Parameters::setBeta(beta);
    Parameters::setFilePath(pwd);
    Parameters::setNSpatial(N);
    Parameters::setNTemporal(N_T);
    Parameters::setNCf(NCf);
    Parameters::setNCor(NCor);
    Parameters::setNTherm(NTherm);
    Parameters::setNUpdates(NUpdates);
    Parameters::setNFlows(NFlows);
    // Setting parallel variables
    Parallel::Communicator::setNumproc(numprocs);
    Parallel::Communicator::setProcessRank(processRank);
    // Main program part
    SU3MatrixGenerator SU3Gen(SU3Eps, seed);
    if (runUnitTestsFlag) {
        if (processRank == 0) runUnitTests(&SU3Gen, runVerboseUnitTestsFlag);
        MPI_Finalize();
        exit(0);
    }
    Plaquette G;
    WilsonGaugeAction S(beta);

    System pureGauge(metropolisSeed, &G, &S);
    if (numberOfArguments > 14) pureGauge.setSubLatticeDimensions(NSub);
    pureGauge.latticeSetup(&SU3Gen, hotStart);
//    pureGauge.setProgramPath(pwd);
//    pureGauge.setConfigBatchName(batchName);
    // ADD VERBOSE ARGUMENT? Or not, more info is always good...?
    pureGauge.printRunInfo(true); // Always print run info, so make optional to turn off
    pureGauge.runMetropolis(storeThermalizationPlaquettes, writeConfigsToFile);
    pureGauge.runBasicStatistics(); // Will always run basic statistics, so turn off!
    pureGauge.printAcceptanceRate(); // Make optional, will always print run info
    pureGauge.writeDataToFile(); // MAKE OPTIONAL ARGUEMENT, that is will always be executed, but can be turned off by calling this function explicitly and setting it to false

    // Finalizing and printing time taken
    programEnd = steady_clock::now();
    programTime = duration_cast<duration<double>>(programEnd - programStart);
    if (processRank == 0) cout << "Program complete. Time used: " << double(programTime.count())/3600.0 << " hours (" << programTime.count() << " seconds)" << endl;

    MPI_Finalize();
    return 0;
}

void runUnitTests(SU3MatrixGenerator *SU3Gen, bool verbose)
{
//    runBoolTest(1e9);
    TestSuite unitTester;
    unitTester.setRNG(SU3Gen);
    unitTester.runFullTestSuite(verbose);
//    SU3BaseTests();
//    runMatrixPerformanceTest(0.24,std::time(nullptr),1e7,true,false);
}
