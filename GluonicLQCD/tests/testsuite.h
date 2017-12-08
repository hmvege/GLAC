#ifndef TESTSUITE_H
#define TESTSUITE_H

#include "math/matrices/su3matrixgenerator.h"
#include "math/latticemath.h"
#include "math/lattice.h"

class TestSuite
{
private:
    // SU3 Test variables
    SU3 U1, U2, U3, UAdd, USub, UMul, UC, UT, UCT, UTrace, U_RST;

    // SU2 Test variables
    SU2 s1, s2, s3, sAdd, sSub, sMul, sC, sT, sCT, s_r, s_s, s_t;

    // Complex test variables
    complex z1, z2, zAdd, zSub, zMul, zDiv, zConj, zSetToMinus;
    double zNorm, zNormSquared;

    // Lattice test variables. Since only the matrix changes, we can reuse previous variables
    Lattice<SU3> latticeSU3;
    Lattice<complex> latticeComplex;
    double latticeDoubleValue;
    Lattice<double> latticeDouble;

    // Limit we demand matrix properties to be retained at
    double m_eps = 2*1e-14;
    double m_tracedMatrix;

    // Verbose storage
    bool m_verbose = false;

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;

    // Basic matrix property testers
    bool operationSU2Test(SU2 results, SU2 solution, std::string operation);
    bool operationSU3Test(SU3 results, SU3 solution, std::string operation);

    // Inline matrix comparing functions
    inline bool compareSU3(SU3 A, SU3 B);
    inline bool compareSU2(SU2 A, SU2 B);
    inline bool compareComplex(complex a, complex b);

    // Inline complex dot product function
    inline complex dot(complex * a, complex * b);
    inline complex dot2(complex * a, complex * b);

    // Basic SU3 matrix operation tests
    bool testSU3Addition();
    bool testSU3Subtraction();
    bool testSU3Multiplication();
    bool testSU3Conjugation();
    bool testSU3Transpose();
    bool testSU3ComplexConjugation();

    // Basic SU2 matrix operation tests
    bool testSU2Addition();
    bool testSU2Subtraction();
    bool testSU2Multiplication();
    bool testSU2Conjugation();
    bool testSU2Transpose();
    bool testSU2ComplexConjugation();

    // SU2 matrix tests
    bool testSU2Hermicity();
    bool checkSU2Hermicity(SU2 H);
    bool testSU2Orthogonality();
    bool checkSU2Orthogonality(SU2 H);
    bool testSU2Norm();
    bool checkSU2Norm(SU2 H);
    bool testSU2Determinant();
    bool checkSU2Determinant(SU2 H);

    // SU3 matrix tests and their sub-functions
    bool testSU3Hermicity();
    bool checkSU3Hermicity(SU3 H);
    bool testSU3Orthogonality();
    bool checkSU3Orthogonality(SU3 H);
    bool testSU3Norm();
    bool checkSU3Norm(SU3 H);
    bool testSU3Determinant();
    bool checkSU3Determinant(SU3 H);

    // Complex class tests
    bool testComplexAddition();
    bool testComplexSubtraction();
    bool testComplexMultiplication();
    bool testComplexDivision();
    bool testComplexConjugation();
    bool testComplexNorm();
    bool testComplexNormSquared();
    bool testComplexSetToMinus();

    // Lattice class tests
    bool testLatticeAddition(); // Test su3, complex and double
    bool testLatticeSubtraction();
    bool testLatticeMultiplication();
    bool testLatticeDivision();
    bool testLatticeRealTrace();
    bool testLatticeImagTrace();
    bool testLatticeSubtractReal();
    bool testLatticeSubtractImag();
    bool testLatticeSum();
    bool testLatticeSumRealTrace();
    bool testLatticeSumRealTraceMultiplication();
    bool testLatticeInverse();


    // Other tests
//    bool testGaugeInvariance();
    bool testSU3TraceMultiplication();
    bool testRSTMultiplication();


public:
    TestSuite();

    void runFullTestSuite(bool verbose);
    bool runComplexTests();
    bool runLatticeTests();
    bool runSU2Tests();
    bool runSU3Tests();
    bool run2x2MatrixTests();
    bool run3x3MatrixTests();
    bool runFunctionsTest();
    void testMatrix(SU3 X);
};

#endif // TESTSUITE_H
