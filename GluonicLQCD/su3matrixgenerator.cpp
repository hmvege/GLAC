#include "su3matrixgenerator.h"
#include <random>
#include <iomanip>
#include <iostream>

#include "su3.h"
#include "complex.h"

// TEMP:
#include "unittests.h"

using std::cout;
using std::endl;

/*
 * Function for generating random SU3 matrices
 */

SU3MatrixGenerator::SU3MatrixGenerator(double eps, std::mt19937_64 &gen, std::uniform_real_distribution<double> &randDistr)
{
    epsilon = eps;
    generator = gen;
    uniform_distribution = randDistr;
}

SU3MatrixGenerator::~SU3MatrixGenerator()
{

}

SU3 SU3MatrixGenerator::generate()
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2
     * 3 4 5
     * 6 7 8
     */
    SU3 H;
    // Generating two random columns
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            H[3*i + j].re = uniform_distribution(generator);
            H[3*i + j].im = uniform_distribution(generator);
        }
    }
    // Using Gram-Schmitt to generate next column
    complex projectionFactor(0,0);
    complex nom(0,0);
    complex denom(0,0);
    for (int i = 0; i < 3; i++)
    {
        nom += H[3*i]*H[3*i+1];
        denom += H[3*i]*H[3*i];
    }
    projectionFactor = nom/denom;
    for (int i = 0; i < 3; i++)
    {
        H[3*i+1] -= projectionFactor*H[3*i];
    }
    // Taking cross product to produce the last column of our matrix
    H[2] = H[3]*H[7] - H[6]*H[4];
    H[5] = H[6]*H[1] - H[0]*H[7];
    H[8] = H[0]*H[4] - H[3]*H[1];

    return H;
}

void SU3MatrixGenerator::generateHermitian()
{

}
