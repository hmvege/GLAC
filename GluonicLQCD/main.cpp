#include <iostream>
#include <ctime>
#include "metropolis.h"
#include "actions/action.h"
#include "su3matrixgenerator.h"

#include "unittests.h"

using std::cout;
using std::endl;

/*
 * TODO:
 * [ ] Change to updating random matrices by X=RST
 * [ ] Add plaquette correlator
 *
 *
 */

int main()
{
    int N           = 8;            // Points for each lattice dimension
    double L        = 2.0;          // Fermi
    int NTherm      = 10;           // Number of times we are to thermalize, that is NTherm * NCor
    int NCor        = 20;           // Only keeping every 20th path
    int NCf         = 1e3;          // Number of configurations to retrieve
    double a        = L/double(N);  // Lattice spacing
    double g        = 5.5;          // Coupling
    double beta     = 6/(g*g);      // Should be
    double SU3Eps   = 0.24;          // Epsilon used for generating SU(3) matrices

    std::mt19937_64 gen(std::time(nullptr));
    std::uniform_real_distribution<double> uni_dist(-SU3Eps,SU3Eps);
    SU3MatrixGenerator SU3Gen(SU3Eps, gen, uni_dist);
    Action S(N,a,beta);
    Metropolis gluon(N, NCf, NCor, NTherm, a, L);
    gluon.setAction(&S);
    gluon.latticeSetup(&SU3Gen);
    gluon.runMetropolis();

    cout << "Program done." << endl;
    return 0;
}
