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
    for (int i = 0; i < 3; i++)
    {
        H[3*i + i].re = 1;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            H[3*i + j].re = uniform_distribution(generator);
            H[3*i + j].im = uniform_distribution(generator);
        }
    }
    // Normalizing first column
    double columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H[3*i].norm();
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[3*i] /= columnLength;
    }
    // Using Gram-Schmitt to generate next column
    complex projectionFactor;
    for (int i = 0; i < 3; i++)
    {
        projectionFactor += H[3*i+1]*H[3*i].c();
    }
    for (int i = 0; i < 3; i++)
    {
        H[3*i+1] = H[3*i+1] - H[3*i]*projectionFactor ;
    }
    // Normalizing second column
    columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H[3*i+1].norm();
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[3*i+1] /= columnLength;
    }
    // Taking cross product to produce the last column of our matrix
    H[2] = H[3].c()*H[7].c() - H[6].c()*H[4].c();
    H[5] = H[6].c()*H[1].c() - H[0].c()*H[7].c();
    H[8] = H[0].c()*H[4].c() - H[3].c()*H[1].c();

//    testOrthogonality(H,true);
//    testNorm(0,H);
//    testNorm(1,H);
//    testNorm(2,H);
//    testHermicity(H,true);
//    exit(1);

    return H;
}


void SU3MatrixGenerator::generateHermitian()
{

}

void SU3MatrixGenerator::GramSchmitt()
{

}
