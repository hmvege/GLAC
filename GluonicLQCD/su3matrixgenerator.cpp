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
    double epsilon = 0.24;
    SU3 H,M;
    for (int i = 0; i < 3; i++)
    {
        H[3*i + i].re = 1;
    }
    // Generating two random columns
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < i+1; j++)
        {
            H[3*i + j].re = uniform_distribution(generator);
            H[3*i + j].im = uniform_distribution(generator);
        }
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < i; j++)
        {
            H[3*j + i].re = H[3*i + j].re;
            H[3*j + i].im = - H[3*i + j].im;
        }
    }
//    for (int i = 0; i < 3; i++)
//    {
//        M[3*i+i].re = 1;
//        for (int j = 0; j < 3; j++)
//        {
//            M[3*i + j].re -= epsilon*H[3*i+j].im;
//            M[3*i + j].im += epsilon*H[3*i+j].re;
//        }
//    }
//    H.copy(M);
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
    H[2] = H[3]*H[7] - H[6]*H[4];
    H[5] = H[6]*H[1] - H[0]*H[7];
    H[8] = H[0]*H[4] - H[3]*H[1];

    testOrthogonality(H,true);
    testNorm(0,H);
    testNorm(1,H);
    testNorm(2,H);
    SU3 HInv;
    HInv.copy(H);
    HInv.transpose();
    HInv.conjugate();
    SU3 I;
    I = H*HInv;

    cout << endl;
    H.print();
    cout << endl;
    HInv.print();
    cout << endl;
    I.print();
    cout << endl;

    exit(1);
    return H;
}

//SU3 SU3MatrixGenerator::generateInverse(SU3 H)
//{
//    SU3 HInv;

//    return HInv;
//}

void SU3MatrixGenerator::generateHermitian()
{

}
