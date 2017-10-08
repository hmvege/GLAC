#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include "su3.h"
#include "su2.h"
#include "complex.h"

class SU3MatrixGenerator
{
private:
    double m_epsilon;
    double m_epsilonSquared;
    double m_sqrtOneMinusEpsSquared;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_distribution;
    std::uniform_real_distribution<double> SU2_uniform_distribution;
    SU3 RSTMatrixMultiplication(SU2 r, SU2 s, SU2 t);
    SU3 RSTMatrixMultiplicationInverse(SU2 r, SU2 s, SU2 t);

    // Used for creating a random matrix
    SU3 H;
    // Used for creating a random matrix close to unity
//    SU3 R, S, T;
    SU3 X;
    SU2 r,s,t;
//    double r0s3, r0s2, r1s3, r1s2, r4s2, r5s2, r4s3, r5s3;
    double rs[8];

    // Used for generating SU2 matrices
    SU2 U;
    double _r[4];
    double _x[4];
    double _rNorm = 0;

    // Pauli matrices
    SU2 sigma[3];
public:
    SU3MatrixGenerator(double eps, double seed);
    ~SU3MatrixGenerator();
    SU3 generateRandom();
    SU3 generateRST();
    SU2 generateSU2();

    // Setters
    void setEpsilon(double epsilon);

    // Getters
    double getEpsilon() { return m_epsilon; }

    // Testers
    SU3 testRSTMultiplication(SU2 r, SU2 s, SU2 t) { return RSTMatrixMultiplication(r,s,t); }
};

#endif // SU3MATRIXGENERATOR_H
