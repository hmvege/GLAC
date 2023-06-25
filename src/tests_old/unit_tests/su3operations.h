#ifndef SU3OPERATIONS_H
#define SU3OPERATIONS_H

#include "testcore.h"

class SU3Operations : public TestCore
{
private:
    // Matrix operation tester
    bool operationSU3Test(SU3 results, SU3 solution, std::string operation);

    // Basic SU3 matrix operation tests
    bool testSU3Addition();
    bool testSU3Subtraction();
    bool testSU3Multiplication();
    bool testSU3Conjugation();
    bool testSU3Transpose();
    bool testSU3ComplexConjugation();

    // SU3 matrix tests and their sub-functions
    bool testSU3Hermicity();
    bool checkSU3Hermicity(SU3 H);
    bool testSU3Orthogonality();
    bool checkSU3Orthogonality(SU3 H);
    bool testSU3Norm();
    bool checkSU3Norm(SU3 H);
    bool testSU3Determinant();
    bool checkSU3Determinant(SU3 H);
    bool testSU3Trace();
    bool testHermitian();
    bool testAntiHermitian();

    // functions.h tests
    bool testSU3TraceMultiplication();

    // Function holding all basic matrix operation tests
    bool run3x3MatrixTests();
public:
    SU3Operations();

    bool runSU3Tests();

    // Checks that the properties of an SU3 matrix is upheld.
    void testMatrix(SU3 X);
};

#endif // SU3OPERATIONS_H
