#include "su3matrixgenerator.h"
#include <random>
#include <iomanip>
#include <iostream>

#include "su3.h"
#include "complex.h"

#include "su2.h"
#include "functions.h"

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
    std::uniform_real_distribution<double> uni_dist_SU2(-0.5,0.5);
    epsilon = eps;
    generator = gen;
    uniform_distribution = uni_dist;
    SU2_uniform_distribution = uni_dist_SU2;
    // Setting up Pauli matrices
    sigma = new SU2[3];
    sigma[0].mat[1].re = 1;
    sigma[0].mat[2].re = 1;
    sigma[1].mat[1].im = -1;
    sigma[1].mat[2].im = 1;
    sigma[2].mat[0].re = 1;
    sigma[2].mat[3].re = -1;
    su2Identity[0].re = 1;
    su2Identity[3].re = 1;
}

SU3MatrixGenerator::~SU3MatrixGenerator()
{
    delete [] sigma;
}

SU3 SU3MatrixGenerator::generateRandom()
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

SU2 SU3MatrixGenerator::generateSU2()
{
    SU2 U;
    double *r = new double[4];
    double *x = new double[4];
    double rNorm = 0;
//    double eps = 0.01;
    // Generating 4 r random numbers
    for (int i = 0; i < 4; i++) {
        r[i] = SU2_uniform_distribution(generator);
    }
    // Populating the x-vector
    if (r[0] < 0) {
        x[0] = - sqrt(1 - epsilon*epsilon);
    }
    else {
        x[0] = sqrt(1 - epsilon*epsilon);
    }
    for (int i = 0; i < 3; i++) {
        rNorm += r[i+1]*r[i+1];
    }
    rNorm = sqrt(rNorm);
    for (int i = 0; i < 3; i++) {
        x[i+1] = epsilon*r[i+1]/rNorm;
    }
    // Imposing unity
    rNorm = 0;
    for (int i = 0; i < 4; i++) {
        rNorm += x[i]*x[i];
    }
    rNorm = sqrt(rNorm);
    for (int i = 0; i < 4; i++) {
        x[i] /= rNorm;
    }
    // Generating the SU2 matrix close to unity
    U += su2Identity*x[0];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            U.mat[j].re += - x[i+1]*sigma[i].mat[j].im;
            U.mat[j].im += x[i+1]*sigma[i].mat[j].re;
        }
    }
    return U;
}

SU3 SU3MatrixGenerator::generateRST()
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
    SU3 X, R, S, T;
    SU2 r,s,t;
    // Generates SU2 matrices
    r = generateSU2();
    s = generateSU2();
    t = generateSU2();
    // Populates R,S,T matrices
    R.mat[0] = r.mat[0];
    R.mat[1] = r.mat[1];
    R.mat[3] = r.mat[2];
    R.mat[4] = r.mat[3];
    R.mat[8].re = 1;
    S.mat[0] = s.mat[0];
    S.mat[2] = s.mat[1];
    S.mat[6] = s.mat[2];
    S.mat[8] = s.mat[3];
    S.mat[4].re = 1;
    T.mat[4] = t.mat[0];
    T.mat[5] = t.mat[1];
    T.mat[7] = t.mat[2];
    T.mat[8] = t.mat[3];
    T.mat[0].re = 1;
    // Creates the return matrix, which is close to unity
    X = R*S*T;
    // Equal probability of returning X and X inverse
    if (SU2_uniform_distribution(generator) < 0) {
        return inverse(X);
    } else {
        return X;
    }
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
