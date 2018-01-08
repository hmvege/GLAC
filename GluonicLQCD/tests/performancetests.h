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

public:
    PerformanceTests();

    void testExponentiation(int NTests);
};

#endif // PERFORMANCETESTS_H
