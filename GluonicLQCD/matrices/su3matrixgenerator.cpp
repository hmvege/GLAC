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
    m_epsilon = eps;
    m_epsilonSquared = eps*eps;
    m_sqrtOneMinusEpsSquared = sqrt(1 - m_epsilonSquared);
    generator = gen;
    uniform_distribution = uni_dist;
    SU2_uniform_distribution = uni_dist_SU2;
    // Setting up Pauli matrices
//    sigma = new SU2[3]; // WOULD THIS BE FASTER IF I USED MEMORY ON STACK?

    sigma[0].mat[2] = 1;
    sigma[0].mat[4] = 1;
    sigma[1].mat[3] = -1;
    sigma[1].mat[5] = 1;
    sigma[2].mat[0] = 1;
    sigma[2].mat[6] = -1;
}

SU3MatrixGenerator::~SU3MatrixGenerator()
{
//    cout << "oops"<<endl;
//    delete [] sigma;
}

void SU3MatrixGenerator::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
    m_epsilonSquared = epsilon;
    m_sqrtOneMinusEpsSquared = sqrt(1 - m_epsilonSquared);
}

SU3 SU3MatrixGenerator::generateRandom()
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
    H.zeros();
    // Populating the matrix
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            H[6*i + 2*j] = uniform_distribution(generator);
            H[6*i + 2*j + 1] = uniform_distribution(generator);
        }
    }
    // Normalizing first column
    double columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H[6*i]*H[6*i] + H[6*i+1]*H[6*i+1];
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[6*i] /= columnLength;
        H[6*i+1] /= columnLength;
    }
    // Using Gram-Schmitt to orthogonalize the second column
    double projectionFactor[2];
    projectionFactor[0] = 0;
    projectionFactor[1] = 0;
    for (int i = 0; i < 3; i++)
    {
        projectionFactor[0] += H[6*i+2]*H[6*i] + H[6*i+3]*H[6*i+1];
        projectionFactor[1] += H[6*i+3]*H[6*i] - H[6*i+2]*H[6*i+1];
    }
    for (int i = 0; i < 3; i++)
    {
        H[6*i+2] -= H[6*i]*projectionFactor[0] - H[6*i+1]*projectionFactor[1];  // Need complex multiplication here
        H[6*i+3] -= H[6*i]*projectionFactor[1] + H[6*i+1]*projectionFactor[0];  // Need complex multiplication here
    }
    // Normalizing second column
    columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H[6*i+2]*H[6*i+2] + H[6*i+3]*H[6*i+3];
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[6*i+2] /= columnLength;
        H[6*i+3] /= columnLength;
    }
    // Taking cross product to produce the last column of our matrix
    H[4] = H[8]*H[12] - H[9]*H[13] - H[14]*H[6] + H[15]*H[7];
    H[5] = H[14]*H[7] + H[15]*H[6] - H[8]*H[13] - H[9]*H[12];
    H[10] = H[14]*H[0] - H[15]*H[1] - H[2]*H[12] + H[3]*H[13];
    H[11] = H[2]*H[13] + H[3]*H[12] - H[14]*H[1] - H[15]*H[0];
    H[16] = H[2]*H[6] - H[3]*H[7] - H[8]*H[0] + H[9]*H[1];
    H[17] = H[8]*H[1] + H[9]*H[0] - H[2]*H[7] - H[3]*H[6];

//    testNorm(0,&H);
//    testNorm(1,&H);
//    testNorm(2,&H);
//    testOrthogonality(H,true);
//    testHermicity(H,false);
//    exit(1);

    return H;
}

SU2 SU3MatrixGenerator::generateSU2()
{
    U.zeros(); // Not needed?
    _rNorm = 0;
    // Generating 4 r random numbers
    for (int i = 0; i < 4; i++) {
        _r[i] = SU2_uniform_distribution(generator);
    }
    // Populating the x-vector
    if (_r[0] < 0) {
        _x[0] = - m_sqrtOneMinusEpsSquared;
    }
    else {
        _x[0] = m_sqrtOneMinusEpsSquared;
    }
    for (int i = 1; i < 4; i++) {
        _rNorm += _r[i]*_r[i];
    }
    _rNorm = sqrt(_rNorm);
    for (int i = 1; i < 4; i++) {
        _x[i] = m_epsilon*_r[i]/_rNorm;
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
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            U[2*j] = U[2*j] - _x[i+1]*sigma[i].mat[2*j+1];
            U[2*j+1] = U[2*j+1] + _x[i+1]*sigma[i].mat[2*j];
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
    // Generates SU2 matrices
    r = generateSU2();
    s = generateSU2();
    t = generateSU2();
    // Populates R,S,T matrices
//    R.zeros(); // No need to populate to unity, as that is covered when populating the SU3 matrices with SU2 matrices.
//    S.zeros();
//    T.zeros();
//    R.mat[0] = r.mat[0];
//    R.mat[1] = r.mat[1];
//    R.mat[2] = r.mat[2];
//    R.mat[3] = r.mat[3];
//    R.mat[6] = r.mat[4];
//    R.mat[7] = r.mat[5];
//    R.mat[8] = r.mat[6];
//    R.mat[9] = r.mat[7];
//    R.mat[16] = 1;
//    S.mat[0] = s.mat[0];
//    S.mat[1] = s.mat[1];
//    S.mat[4] = s.mat[2];
//    S.mat[5] = s.mat[3];
//    S.mat[12] = s.mat[4];
//    S.mat[13] = s.mat[5];
//    S.mat[16] = s.mat[6];
//    S.mat[17] = s.mat[7];
//    S.mat[8] = 1;
//    T.mat[8] = t.mat[0];
//    T.mat[9] = t.mat[1];
//    T.mat[10] = t.mat[2];
//    T.mat[11] = t.mat[3];
//    T.mat[14] = t.mat[4];
//    T.mat[15] = t.mat[5];
//    T.mat[16] = t.mat[6];
//    T.mat[17] = t.mat[7];
//    T.mat[0] = 1;
//    X.zeros();
    // CLEAN UP THIS! ALSO; MAKE DIRECT INVERSE METHOD!!
//    X[0] = - r[1]*s[1] + r[0]*s[0];
//    X[1] = + r[1]*s[0] + r[0]*s[1];
//    X[2] = - r[1]*s[3]*t[4] - r[1]*s[2]*t[5] - r[0]*s[3]*t[5] + r[0]*s[2]*t[4] - r[3]*t[1] + r[2]*t[0];
//    X[3] = - r[1]*s[3]*t[5] + r[1]*s[2]*t[4] + r[0]*s[3]*t[4] + r[0]*s[2]*t[5] + r[3]*t[0] + r[2]*t[1];
//    X[4] = - r[1]*s[3]*t[6] - r[1]*s[2]*t[7] - r[0]*s[3]*t[7] + r[0]*s[2]*t[6] - r[3]*t[3] + r[2]*t[2];
//    X[5] = - r[1]*s[3]*t[7] + r[1]*s[2]*t[6] + r[0]*s[3]*t[6] + r[0]*s[2]*t[7] + r[3]*t[2] + r[2]*t[3];
//    X[6] = - r[5]*s[1] + r[4]*s[0];
//    X[7] = + r[5]*s[0] + r[4]*s[1];
//    X[8] = - r[5]*s[3]*t[4] - r[5]*s[2]*t[5] - r[4]*s[3]*t[5] + r[4]*s[2]*t[4] - r[7]*t[1] + r[6]*t[0];
//    X[9] = - r[5]*s[3]*t[5] + r[5]*s[2]*t[4] + r[4]*s[3]*t[4] + r[4]*s[2]*t[5] + r[7]*t[0] + r[6]*t[1];
//    X[10] = - r[5]*s[3]*t[6] - r[5]*s[2]*t[7] - r[4]*s[3]*t[7] + r[4]*s[2]*t[6] - r[7]*t[3] + r[6]*t[2];
//    X[11] = - r[5]*s[3]*t[7] + r[5]*s[2]*t[6] + r[4]*s[3]*t[6] + r[4]*s[2]*t[7] + r[7]*t[2] + r[6]*t[3];
//    X[12] = + s[4];
//    X[13] =  s[5];
//    X[14] = - s[7]*t[5] + s[6]*t[4];
//    X[15] = + s[7]*t[4] + s[6]*t[5];
//    X[16] = - s[7]*t[7] + s[6]*t[6];
//    X[17] = + s[7]*t[6] + s[6]*t[7];

//    SU3 H = R*S*T;
//    testOrthogonality(H,false);
//    testNorm(0,H);
//    testNorm(1,H);
//    testNorm(2,H);
//    testHermicity(H,false);
//    cout << "Completed SU3 RST test" << endl;
//    exit(1);

    if (SU2_uniform_distribution(generator) < 0) {
        return RSTMatrixMultiplication(r,s,t).inv();
//        return (R*S*T).inv(); // THIS CAN BE MADE INTO TWO RST FUNCTIONS WITH 20 CFLOPS!!
//        return X.inv(); // THIS CAN BE MADE INTO TWO RST FUNCTIONS WITH 20 CFLOPS!!
    } else {
        return RSTMatrixMultiplication(r,s,t);
//        return X;
//        return R*S*T;
    }
}

SU3 SU3MatrixGenerator::RSTMatrixMultiplication(SU2 r, SU2 s, SU2 t)
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
//    H.mat[]
    X.zeros();
    X[0] = r[0]*s[0] - r[1]*s[1];
    X[1] = r[1]*s[0] + r[0]*s[1];
    X[2] = r[0]*s[2]*t[4] - r[3]*t[1] + r[2]*t[0] - r[1]*s[3]*t[4] - r[1]*s[2]*t[5] - r[0]*s[3]*t[5];
    X[3] = r[1]*s[2]*t[4] + r[0]*s[3]*t[4] + r[0]*s[2]*t[5] + r[3]*t[0] + r[2]*t[1] - r[1]*s[3]*t[5];
    X[4] = r[0]*s[2]*t[6] - r[3]*t[3] + r[2]*t[2] - r[1]*s[3]*t[6] - r[1]*s[2]*t[7] - r[0]*s[3]*t[7];
    X[5] = r[1]*s[2]*t[6] + r[0]*s[3]*t[6] + r[0]*s[2]*t[7] + r[3]*t[2] + r[2]*t[3] - r[1]*s[3]*t[7];
    X[6] = r[4]*s[0] - r[5]*s[1];
    X[7] = r[5]*s[0] + r[4]*s[1];
    X[8] = r[4]*s[2]*t[4] - r[7]*t[1] + r[6]*t[0] - r[5]*s[3]*t[4] - r[5]*s[2]*t[5] - r[4]*s[3]*t[5];
    X[9] = r[5]*s[2]*t[4] + r[4]*s[3]*t[4] + r[4]*s[2]*t[5] + r[7]*t[0] + r[6]*t[1] - r[5]*s[3]*t[5];
    X[10] = r[4]*s[2]*t[6] - r[7]*t[3] + r[6]*t[2] - r[5]*s[3]*t[6] - r[5]*s[2]*t[7] - r[4]*s[3]*t[7];
    X[11] = r[5]*s[2]*t[6] + r[4]*s[3]*t[6] + r[4]*s[2]*t[7] + r[7]*t[2] + r[6]*t[3] - r[5]*s[3]*t[7];
    X[12] = s[4];
    X[13] = s[5];
    X[14] = s[6]*t[4] - s[7]*t[5];
    X[15] = s[7]*t[4] + s[6]*t[5];
    X[16] = s[6]*t[6] - s[7]*t[7];
    X[17] = s[7]*t[6] + s[6]*t[7];
    return X;
}


SU3 SU3MatrixGenerator::RSTMatrixMultiplicationInverse(SU2 r, SU2 s, SU2 t)
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
//    H.mat[]
}
