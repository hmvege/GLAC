#include <iostream>
#include <ctime>
#include <chrono>
#include "system.h"
#include "actions/action.h"
#include "actions/wilsongaugeaction.h"
#include "correlators/plaquette.h"
#include "matrices/su3matrixgenerator.h"

#include "unittests.h"

using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

/*
 * TODO:
 * [ ] Switch to CORRECT method syntax, foo --> m_foo
 * [ ] Check that the lattice is gauge invariant: M^-1 * U * M, see Gattinger intro on how to make gauge fields gauge invariant!
 */

int main()
{
    int N           = 8;            // Points for each lattice dimension, twice the points in time dimension
    int N_T         = 8;            // Time dimension
    int NTherm      = 100;           // Number of times we are to thermalize. Default 200-300.
    int NCor        = 10;           // Only keeping every 20th path
    int NCf         = 5;          // Number of configurations to retrieve
    double beta     = 6;            // Should be
    double SU3Eps   = 0.24;         // Epsilon used for generating SU(3) matrices
    double seed     = std::time(nullptr);
    double metropolisSeed = std::time(nullptr) + 1;
    bool storeThermalizationPlaquettes = false;
    std::string filename = "TEST_SCALAR_RUN1";

//    testMatrixSU3Properties();
//    runMatrixPerformanceTest(SU3Eps, seed, 1e8,false,true);
//    exit(1);
    // Timers

    steady_clock::time_point programStart, programEnd;
    duration<double> programTime;
    programStart = steady_clock::now();


    SU3MatrixGenerator SU3Gen(SU3Eps, seed);
    Plaquette G(N, N_T);
    WilsonGaugeAction S(N, N_T, beta);
    System pureGauge(N, N_T, NCf, NCor, NTherm, metropolisSeed, &G, &S);
    pureGauge.latticeSetup(&SU3Gen);
    pureGauge.runMetropolis(storeThermalizationPlaquettes);
    pureGauge.getStatistics();
    pureGauge.printAcceptanceRate();
//    pureGauge.writeConfigurationToFile(filename); // Only needed when generating configs to parallel processors
    pureGauge.writeDataToFile(filename);

    programEnd = steady_clock::now();
    programTime = duration_cast<duration<double>>(programEnd - programStart);
    cout << "Program complete. Time used: " << programTime.count() << endl;
    return 0;
}
