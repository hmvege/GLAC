#ifndef TESTSUITE_H
#define TESTSUITE_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "math/matrices/su2.h"
#include "math/matrices/su3.h"
#include "math/matrices/su3matrixgenerator.h"
#include "observables/plaquette.h"
#include "math/links.h"
#include "math/complex.h"
#include "math/functions.h"

class TestSuite
{
private:
    SU3 U1, U2, U3, UAdd, USub, UMul, UC, UT, UCT, UTrace, U_RST;
    SU2 s1, s2, s3, sAdd, sSub, sMul, sC, sT, sCT, s_r, s_s, s_t;

    double m_eps = 2*1e-14;
    double m_tracedMatrix;

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;

    // Basic matrix property testers
    bool operationSU2Test(bool verbose, SU2 results, SU2 solution, std::string operation);
    bool operationSU3Test(bool verbose, SU3 results, SU3 solution, std::string operation);

    // Inline matrix comparing functions
    inline bool compareSU3(SU3 A, SU3 B);
    inline bool compareSU2(SU2 A, SU2 B);

    // Inline complex dot product function
    inline complex dot(complex * a, complex * b);
    inline complex dot2(complex * a, complex * b);

    // Basic SU3 matrix operation tests
    bool testSU3Addition(bool verbose);
    bool testSU3Subtraction(bool verbose);
    bool testSU3Multiplication(bool verbose);
    bool testSU3Conjugation(bool verbose);
    bool testSU3Transpose(bool verbose);
    bool testSU3ComplexConjugation(bool verbose);

    // Basic SU2 matrix operation tests
    bool testSU2Addition(bool verbose);
    bool testSU2Subtraction(bool verbose);
    bool testSU2Multiplication(bool verbose);
    bool testSU2Conjugation(bool verbose);
    bool testSU2Transpose(bool verbose);
    bool testSU2ComplexConjugation(bool verbose);

    // SU2 matrix tests
    bool testSU2Hermicity(bool verbose);
    bool checkSU2Hermicity(bool verbose, SU2 H);
    bool testSU2Orthogonality(bool verbose);
    bool checkSU2Orthogonality(bool verbose, SU2 H);
    bool testSU2Norm(bool verbose);
    bool checkSU2Norm(bool verbose, SU2 H);
    bool testSU2Determinant(bool verbose);
    bool checkSU2Determinant(bool verbose, SU2 H);

    // SU3 matrix tests and their sub-functions
    bool testSU3Hermicity(bool verbose);
    bool checkSU3Hermicity(bool verbose, SU3 H);
    bool testSU3Orthogonality(bool verbose);
    bool checkSU3Orthogonality(bool verbose, SU3 H);
    bool testSU3Norm(bool verbose);
    bool checkSU3Norm(bool verbose, SU3 H);
    bool testSU3Determinant(bool verbose);
    bool checkSU3Determinant(bool verbose, SU3 H);

    // Other tests
//    bool testGaugeInvariance(bool verbose);
    bool testSU3TraceMultiplication(bool verbose);
    bool testRSTMultiplication(bool verbose);

    // Add complex class tests?

    // Add performance tests?
public:
    TestSuite();

    void runFullTestSuite(bool verbose);
    bool runSU2Tests(bool verbose);
    bool runSU3Tests(bool verbose);
    bool run2x2MatrixTests(bool verbose);
    bool run3x3MatrixTests(bool verbose);
    bool runFunctionsTest(bool verbose);
    void testMatrix(SU3 X, bool verbose);

    // Setters
    void setRNG(SU3MatrixGenerator *SU3Gen);
};

#endif // TESTSUITE_H
