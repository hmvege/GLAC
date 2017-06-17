#include <iostream>
#include <ctime>
#include "metropolis.h"
#include "actions/action.h"
#include "su3matrixgenerator.h"

#include "unittests.h"

using namespace std;

int main()
{
    int N           = 8;            // Points for each lattice dimension
    double L        = 2.0;          // Fermi
    int NTherm      = 10;           // Number of times we are to thermalize, that is NTherm * NCor
    int NCor        = 20;           // Only keeping every 20th path
    int NCf         = 1e5;          // Number of random path or path configurations
    double a        = L/double(N);  // Lattice spacing
    double SU3Eps   = 0.1;          // Epsilon used for generating SU(3) matrices

    std::mt19937_64 gen(std::time(nullptr));
    std::uniform_real_distribution<double> uni_dist(-1,1);
    SU3MatrixGenerator SU3Gen(SU3Eps, gen, uni_dist);
    Metropolis gluon(N, NCf, NCor, NTherm, a, L);
    gluon.latticeSetup(&SU3Gen);

    cout << "Program done." << endl;
    return 0;
}
