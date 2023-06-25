/*!
 * \brief
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef PERFORMANCETESTS_H
#define PERFORMANCETESTS_H

#include "math/matrices/su3matrixgenerator.h"
#include "tests_old/unit_tests/testcore.h"

class PerformanceTests : public TestCore
{
private:
    // Degree of the Taylor polynomial to use in exponentiation
    unsigned int m_NTaylorDegree = 8;

    // Tests for the SU3 exponentiation
    void testExponentiationTime(unsigned int NTests);
    void testExponentiationAccuracy();

    // Tests of the RNGs
    void testRandomGenerators(unsigned int NTests);

    // Tests for the SU3 derivative
    void testDerivativeTimeAndAccuracy(unsigned int NTests);

    // Tests for the SU3 matrix multiplication
    void testMatrixMultiplication();

    // Tests for flow
    void testFlow();

    // Tests for shift
    void testShift();
    void testDoubleShift();
public:
    PerformanceTests();
    ~PerformanceTests();

    void run();
};

#endif // PERFORMANCETESTS_H
