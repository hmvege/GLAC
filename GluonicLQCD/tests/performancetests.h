#ifndef PERFORMANCETESTS_H
#define PERFORMANCETESTS_H

#include "math/matrices/su3matrixgenerator.h"

class PerformanceTests
{
private:
    // SU3 generator
    SU3MatrixGenerator * m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;

    // Degree of the Taylor polynomial to use in exponentiation
    int m_NTaylorDegree = 8;

    // Tests for the SU3 exponentiation
    void testExponentiationTime(int NTests);
    void testExponentiationAccuracy();

    // Tests of the RNGs
    void testRandomGenerators(int NTests);

    // Tests for the SU3 derivative
    void testDerivativeTimeAndAccuracy(int NTests);
public:
    PerformanceTests();

    void run(int NExponentiationTests, int NRandomTests, int NDerivativeTests);
};

#endif // PERFORMANCETESTS_H
