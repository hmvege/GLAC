#include <iostream>
#include <ctime>
#include "system.h"
#include "actions/action.h"
#include "actions/wilsongaugeaction.h"
#include "correlators/plaquette.h"
#include "matrices/su3matrixgenerator.h"

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
 * [ ] Switch to CORRECT method syntax, foo --> m_foo
 * [ ] Check that the lattice is gauge invariant: M^-1 * U * M, see Gattinger intro on how to make gauge fields gauge invariant!
 */

int main()
{
    int N           = 16;            // Points for each lattice dimension, twice the points in time dimension
    int N_T         = 16;            // Time dimension
    int NTherm      = 25;           // Number of times we are to thermalize. Default 200-300.
    int NCor        = 10;           // Only keeping every 20th path
    int NCf         = 5;          // Number of configurations to retrieve
    double beta     = 6;            // Should be
    double SU3Eps   = 0.24;         // Epsilon used for generating SU(3) matrices
    double seed     = std::time(nullptr);
    double metropolisSeed = std::time(nullptr) + 1;
    bool storeThermalizationPlaquettes = true;
    std::string filename = "scalar16cubed16run1";

//    testMatrixSU3Properties();
//    runMatrixPerformanceTest(SU3Eps, seed, 1e8,false,true);
//    exit(1);
    clock_t programStart, programEnd;
    programStart = clock();
    SU3MatrixGenerator SU3Gen(SU3Eps, seed);
    Plaquette G(N, N_T);
    WilsonGaugeAction S(N, N_T, beta);
    System pureGauge(N, N_T, NCf, NCor, NTherm, metropolisSeed, &G, &S);
    pureGauge.latticeSetup(&SU3Gen);
    pureGauge.runMetropolis(storeThermalizationPlaquettes);
    pureGauge.getStatistics();
    pureGauge.printAcceptanceRate();
    pureGauge.writeConfigurationToFile(filename);
    pureGauge.writeDataToFile(filename);

    programEnd = clock();
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}
