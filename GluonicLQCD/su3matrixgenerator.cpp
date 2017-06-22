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

SU3MatrixGenerator::SU3MatrixGenerator(double eps, double seed)
{
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> uni_dist(-eps,eps);
    epsilon = eps;
    generator = gen;
    uniform_distribution = uni_dist;
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
//    for (int i = 0; i < 3; i++)
//    {
//        H[3*i + i].re = 1;
//    }
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
//    testHermicity(H,false);
//    exit(1);

    return H;
}

SU3 SU3MatrixGenerator::generateIdentity()
{
    SU3 H;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i==j)
            {
                H[3*i+j].re = 1;
                H[3*i+j].im = 0;
            }
            else
            {
                H[3*i+j].re = 0;
                H[3*i+j].im = 0;
            }
        }
    }
    return H;
}

SU3 SU3MatrixGenerator::updateMatrix()
{
    /*
     * Generatores a random SU3 matrix close to unity.
     * Index map:
     * r,s,t =
     * 0 1
     * 2 3
     * R,S,T =
     * 0 1 2
     * 3 4 5
     * 6 7 8
     */
    SU3 R, S, T;
    complex * r = new complex[4];
    complex * s = new complex[4];
    complex * t = new complex[4];
    // Generates 3 SU2 matrices
    // Gets a
    r[0].re = uniform_distribution(generator);
    r[0].im = uniform_distribution(generator);
    s[0].re = uniform_distribution(generator);
    s[0].im = uniform_distribution(generator);
    t[0].re = uniform_distribution(generator);
    t[0].im = uniform_distribution(generator);
    // Gets b
    r[1].re = uniform_distribution(generator);
    r[1].im = uniform_distribution(generator);
    r[2] -= r[1].c();
    r[3] = r[0].c();
    s[1].re = uniform_distribution(generator);
    s[1].im = uniform_distribution(generator);
    s[2] -= s[1].c();
    s[3] = s[0].c();
    t[1].re = uniform_distribution(generator);
    t[1].im = uniform_distribution(generator);
    t[2] -= t[1].c();
    t[3] = t[0].c();
    // Insert into R,S and T
    R.mat[0] = r[0];
    R.mat[1] = r[1];
    R.mat[3] = r[2];
    R.mat[4] = r[3];
    R.mat[8].re = 1;
    S.mat[0] = s[0];
    S.mat[2] = s[1];
    S.mat[6] = s[2];
    S.mat[8] = s[3];
    S.mat[4].re = 1;
    T.mat[4] = t[0];
    T.mat[5] = t[1];
    T.mat[7] = t[2];
    T.mat[8] = t[3];
    T.mat[0].re = 1;
    // Returns product
    delete [] r;
    delete [] s;
    delete [] t;
    return R*S*T;
}

void SU3MatrixGenerator::generateHermitian()
{
    cout << "SU3MatrixGenerator::generateHermitian not implemented" << endl;
    exit(1);
}

void SU3MatrixGenerator::GramSchmitt()
{
    cout << "SU3MatrixGenerator::GramSchmitt not implemented" << endl;
    exit(1);
}
