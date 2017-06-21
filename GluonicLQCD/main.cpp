#include <iostream>
#include <ctime>
#include "metropolis.h"
#include "actions/action.h"
#include "actions/wilsongaugeaction.h"
#include "correlators/plaquette.h"
#include "su3matrixgenerator.h"

#include "unittests.h"

using std::cout;
using std::endl;

/*
 * TODO:
 * [x] Add plaquette correlator
 * [x] Make actions more general!! Aka, create a Wilson action
 * [ ] Fix bug that makes me get negative correlators
 * [ ] Change to updating random matrices by X=RST
 * [ ] Switch to method syntax, foo --> m_foo
 * [ ] Create method for saving lattice configuration
 */

int main()
{
    int N           = 4;            // Points for each lattice dimension
    double L        = 2.0;          // Length of lattice in fermi
    int NTherm      = 10;           // Number of times we are to thermalize, that is NTherm * NCor
    int NCor        = 10;           // Only keeping every 20th path
    int NCf         = 1000;          // Number of configurations to retrieve
    double a        = L/double(N);  // Lattice spacing
//    double g        = 5.5;          // Coupling
    double beta     = 6;            // Should be
    double SU3Eps   = 0.24;         // Epsilon used for generating SU(3) matrices

    clock_t programStart, programEnd;
    programStart = clock();
    std::mt19937_64 gen(std::time(nullptr));
    std::uniform_real_distribution<double> uni_dist(-SU3Eps,SU3Eps);
    SU3MatrixGenerator SU3Gen(SU3Eps, gen, uni_dist);
    Plaquette G(N);
    WilsonGaugeAction S(N,beta);
    Metropolis gluon(N, NCf, NCor, NTherm, a, L, &G, &S);
    gluon.latticeSetup(&SU3Gen);
    gluon.runMetropolis();
    gluon.getStatistics();
    gluon.writeDataToFile("../output/gluon_data.txt");

    programEnd = clock();
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}
