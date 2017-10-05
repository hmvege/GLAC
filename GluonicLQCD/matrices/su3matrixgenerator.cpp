#include "su3matrixgenerator.h"
#include <random>
#include <iomanip>
#include <iostream>
#include "functions.h"

// For testing
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
    epsilonSquared = eps*eps;
    generator = gen;
    uniform_distribution = uni_dist;
    SU2_uniform_distribution = uni_dist_SU2;
    // Setting up Pauli matrices
    sigma = new SU2[3]; // WOULD THIS BE FASTER IF I USED MEMORY ON STACK?

    sigma[0].mat[2] = 1;
    sigma[0].mat[4] = 1;
    sigma[1].mat[3] = -1;
    sigma[1].mat[5] = 1;
    sigma[2].mat[0] = 1;
    sigma[2].mat[6] = -1;
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
    cout << "NOT USING REGULAR RANDOM -> SHOULD NOT BE SEEN!" << endl;
    exit(1);
    H.identity();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            H[3*i + j] = uniform_distribution(generator);
            H[3*i + j + 1] = uniform_distribution(generator);
        }
    }
    // Normalizing first column
    double columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
//        columnLength += H[3*i].normSquared();
        columnLength += H.normSquared(6*i);
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[6*i] /= columnLength;
    }
    // Using Gram-Schmitt to generate next column
    complex projectionFactor;
    for (int i = 0; i < 3; i++)
    {
        projectionFactor += complexMultiply(H, H.conjugate(), 6*i+2, 6*i);  // Need complex multiplication here
//        projectionFactor += H[3*i+1]*H[3*i].c();  // Need complex multiplication here
    }
    for (int i = 0; i < 3; i++)
    {
        H[6*i+2] = H[6*i+2] - H[6*i]*projectionFactor.re();  // Need complex multiplication here
        H[6*i+3] = H[6*i+3] - H[6*i+1]*projectionFactor.im();  // Need complex multiplication here
//        H[3*i+1] = H[3*i+1] - H[3*i]*projectionFactor;  // Need complex multiplication here
    }
    // Normalizing second column
    columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H.normSquared(6*i+2);
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[6*i+2] /= columnLength;
    }
    // Taking cross product to produce the last column of our matrix
//    H[2] = H[3].c()*H[7].c() - H[6].c()*H[4].c();  // Need complex multiplication here
//    H[5] = H[6].c()*H[1].c() - H[0].c()*H[7].c();  // Need complex multiplication here
//    H[8] = H[0].c()*H[4].c() - H[3].c()*H[1].c();  // Need complex multiplication here

    H[4] = H[6]*H[14] - H[7]*H[15] + H[12]*H[8] - H[13]*H[9];
    H[5] = -H[6]*H[15] - H[7]*H[14] - H[12]*H[9] - H[13]*H[8];

    H[10] = H[12]*H[2] - H[13]*H[3] + H[0]*H[14] - H[1]*H[15];
    H[11] = - H[12]*H[3] - H[13]*H[2] - H[0]*H[15] - H[1]*H[14];

    H[16] = H[0]*H[8] - H[1]*H[9] + H[6]*H[2] - H[7]*H[3];
    H[17] = - H[0]*H[9] - H[1]*H[8] - H[6]*H[3] - H[7]*H[2];

//    testOrthogonality(H,false);
//    testNorm(0,H);
//    testNorm(1,H);
//    testNorm(2,H);
//    testHermicity(H,false);
//    exit(1);

    return H;
}

SU2 SU3MatrixGenerator::generateSU2()
{
//    SU2 U;
//    double _r[4];
//    double _x[4];
//    double _rNorm = 0;
    U.zeros(); // Not needed?
    _rNorm = 0;
    // Generating 4 r random numbers
    for (int i = 0; i < 4; i++) {
        _r[i] = SU2_uniform_distribution(generator);
    }
    // Populating the x-vector
    if (_r[0] < 0) {
        _x[0] = - sqrt(1 - epsilonSquared);
    }
    else {
        _x[0] = sqrt(1 - epsilonSquared);
    }
//    for (int i = 0; i < 3; i++) {
//        _rNorm += _r[i+1]*_r[i+1];
//    }
    for (int i = 1; i < 4; i++) {
        _rNorm += _r[i]*_r[i];
    }
    _rNorm = sqrt(_rNorm);
//    for (int i = 0; i < 3; i++) {
//        _x[i+1] = epsilon*_r[i+1]/_rNorm;
//    }
    for (int i = 1; i < 4; i++) {
        _x[i] = epsilon*_r[i]/_rNorm;
    }
    // Imposing unity
    _rNorm = 0;
    for (int i = 0; i < 4; i++) {
        _rNorm += _x[i]*_x[i];
    }
    _rNorm = sqrt(_rNorm);
    for (int i = 0; i < 4; i++) {
        _x[i] /= _rNorm;
    }
    // Generating the SU2 matrix close to unity
    U.mat[0] = _x[0]; // same as 1*x0
    U.mat[6] = _x[0]; // same as 1*x0
//    U *= _x[0];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            U[2*j] = U[2*j] - _x[i+1]*sigma[i].mat[2*j+1];
            U[2*j+1] = U[2*j+1] + _x[i+1]*sigma[i].mat[2*j];
        }
    }
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 4; j++) {
//            U.mat[j].setRe(U.mat[j].re() - _x[i+1]*sigma[i].mat[j].im());
//            U.mat[j].setIm(U.mat[j].im() + _x[i+1]*sigma[i].mat[j].re());
//        }
//    }
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
//    SU3 X, R, S, T; // Change to only have 3 matrices instead of 4, such that one is just R=R*S*T?
//    SU2 r,s,t;
    // Generates SU2 matrices
    r = generateSU2();
    s = generateSU2();
    t = generateSU2();
    // Populates R,S,T matrices
    R.zeros(); // No need to populate to unity, as that is covered when populating the SU3 matrices with SU2 matrices.
    S.zeros();
    T.zeros();
    R.mat[0] = r.mat[0];
    R.mat[2] = r.mat[2];
    R.mat[6] = r.mat[4];
    R.mat[8] = r.mat[6];
    R.mat[16] = 1;
    S.mat[0] = s.mat[0];
    S.mat[4] = s.mat[2];
    S.mat[12] = s.mat[4];
    S.mat[16] = s.mat[6];
    S.mat[8] = 1;
    T.mat[8] = t.mat[0];
    T.mat[10] = t.mat[2];
    T.mat[14] = t.mat[4];
    T.mat[16] = t.mat[6];
    T.mat[0] = 1;
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++) {
//            R.mat[3*i + j] = r.mat[2*i + j];
//            S.mat[6*i + 2*j] = s.mat[2*i + j];
//            T.mat[4 + i*3 + j] = t.mat[2*i + j];
//        }
//    }
    // Creates the return matrix, which is close to unity
//    X = R*S*T;
    // Equal probability of returning X and X inverse
    if (SU2_uniform_distribution(generator) < 0) {
        return (R*S*T).inv();
    } else {
        return R*S*T;
    }
}
