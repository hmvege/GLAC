#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include "su3.h"
#include "su2.h"

class SU3MatrixGenerator
{
private:
    // RNG
    std::mt19937_64 m_generator;

    // Random matrix distribution
    std::uniform_real_distribution<double> m_uniform_distribution;

    // RST related distribution
    std::uniform_real_distribution<double> m_SU2_uniform_distribution;
    double m_epsilon;
    double m_epsilonSquared;
    double m_sqrtOneMinusEpsSquared;
    inline SU3 RSTMatrixMultiplication(SU2 r, SU2 s, SU2 t);
    inline SU3 RSTMatrixMultiplicationInverse(SU2 r, SU2 s, SU2 t);

    // Used for creating a random matrix
    SU3 H;

    // Used for creating a random matrix close to unity
    SU3 X;
    SU2 r,s,t;
    double rs[8];

    // Used for generating SU2 matrices
    SU2 U;
    double _r[4];
    double _x[4];
    double _rNorm = 0;

    // Pauli matrices
    SU2 sigma[3];
public:
    SU3MatrixGenerator();
    ~SU3MatrixGenerator();
    SU3 generateRandom();
    SU3 generateRST();
    inline SU2 generateSU2();

    // Testers
    SU3 testRSTMultiplication(SU2 r, SU2 s, SU2 t) { return RSTMatrixMultiplication(r,s,t); }
    SU3 testRSTMultiplicationInverse(SU2 r, SU2 s, SU2 t) { return RSTMatrixMultiplicationInverse(r,s,t); }
};

inline SU2 SU3MatrixGenerator::generateSU2()
{
    U.zeros();
    _rNorm = 0;
    // Generating 4 r random numbers
    for (int i = 0; i < 4; i++) {
        _r[i] = m_SU2_uniform_distribution(m_generator);
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
    _rNorm = sqrt(_rNorm)/m_epsilon;
    for (int i = 1; i < 4; i++) {
        _x[i] = _r[i]/_rNorm;
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

inline SU3 SU3MatrixGenerator::RSTMatrixMultiplication(SU2 r, SU2 s, SU2 t)
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
    // Block one shortenings
    rs[0] = r[0]*s[2] - r[1]*s[3];
    rs[1] = r[1]*s[2] + r[0]*s[3];
    // Block two shortenings
    rs[2] = r[4]*s[2] - r[5]*s[3];
    rs[3] = r[5]*s[2] + r[4]*s[3];
    // Compact RST multiplication
    X[0] = r[0]*s[0] - r[1]*s[1];
    X[1] = r[1]*s[0] + r[0]*s[1];
    X[2] = rs[0]*t[4] - r[3]*t[1] + r[2]*t[0] - rs[1]*t[5];
    X[3] = rs[1]*t[4] + r[3]*t[0] + r[2]*t[1] + rs[0]*t[5];
    X[4] = rs[0]*t[6] - r[3]*t[3] + r[2]*t[2] - rs[1]*t[7];
    X[5] = rs[1]*t[6] + r[3]*t[2] + r[2]*t[3] + rs[0]*t[7];
    X[6] = r[4]*s[0] - r[5]*s[1];
    X[7] = r[5]*s[0] + r[4]*s[1];
    X[8] = rs[2]*t[4] - r[7]*t[1] + r[6]*t[0] - rs[3]*t[5];
    X[9] = rs[3]*t[4] + r[7]*t[0] + r[6]*t[1] + rs[2]*t[5];
    X[10] = rs[2]*t[6] - r[7]*t[3] + r[6]*t[2] - rs[3]*t[7];
    X[11] = rs[3]*t[6] + r[7]*t[2] + r[6]*t[3] + rs[2]*t[7];
    X[12] = s[4];
    X[13] = s[5];
    X[14] = s[6]*t[4] - s[7]*t[5];
    X[15] = s[7]*t[4] + s[6]*t[5];
    X[16] = s[6]*t[6] - s[7]*t[7];
    X[17] = s[7]*t[6] + s[6]*t[7];

    return X;
}


inline SU3 SU3MatrixGenerator::RSTMatrixMultiplicationInverse(SU2 r, SU2 s, SU2 t)
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
    // Block one shortenings
    rs[0] = r[0]*s[2] - r[1]*s[3];
    rs[1] = r[1]*s[2] + r[0]*s[3];
    // Block two shortenings
    rs[2] = r[4]*s[2] - r[5]*s[3];
    rs[3] = r[5]*s[2] + r[4]*s[3];

    // Upper triangular
    X[6]  =   rs[0]*t[4] - r[3]*t[1] + r[2]*t[0] - rs[1]*t[5];
    X[7]  = - rs[1]*t[4] - r[3]*t[0] - r[2]*t[1] - rs[0]*t[5];
    X[12] =   rs[0]*t[6] - r[3]*t[3] + r[2]*t[2] - rs[1]*t[7];
    X[13] = - rs[1]*t[6] - r[3]*t[2] - r[2]*t[3] - rs[0]*t[7];
    X[14] =   rs[2]*t[6] - r[7]*t[3] + r[6]*t[2] - rs[3]*t[7];
    X[15] = - rs[3]*t[6] - r[7]*t[2] - r[6]*t[3] - rs[2]*t[7];
    // Upper triangular
    X[2]  =   r[4]*s[0] - r[5]*s[1];
    X[3]  = - r[5]*s[0] - r[4]*s[1];
    X[4]  =   s[4];
    X[5]  = - s[5];
    X[10] =   s[6]*t[4] - s[7]*t[5];
    X[11] = - s[7]*t[4] - s[6]*t[5];
    // Diagonals
    X[0]  =   r[0]*s[0] - r[1]*s[1];
    X[1]  = - r[1]*s[0] - r[0]*s[1];
    X[8]  =   rs[2]*t[4] - r[7]*t[1] + r[6]*t[0] - rs[3]*t[5];;
    X[9]  = - rs[3]*t[4] - r[7]*t[0] - r[6]*t[1] - rs[2]*t[5];
    X[16] =   s[6]*t[6] - s[7]*t[7];
    X[17] = - s[7]*t[6] - s[6]*t[7];

    return X;
}


#endif // SU3MATRIXGENERATOR_H
