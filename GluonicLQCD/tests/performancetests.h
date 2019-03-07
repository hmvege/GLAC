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
    unsigned int m_NTaylorDegree = 8;

    // Tests for the SU3 exponentiation
    void testExponentiationTime(unsigned int NTests);
    void testExponentiationAccuracy();

    // Tests of the RNGs
    void testRandomGenerators(unsigned int NTests);

    // Tests for the SU3 derivative
    void testDerivativeTimeAndAccuracy(unsigned int NTests);
public:
    PerformanceTests();
    ~PerformanceTests();

    void run();
};

#endif // PERFORMANCETESTS_H
