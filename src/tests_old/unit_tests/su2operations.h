#ifndef SU2OPERATIONS_H
#define SU2OPERATIONS_H

#include "testcore.h"

class SU2Operations : public TestCore
{
private:
    bool operationSU2Test(SU2 results, SU2 solution, std::string operation);

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

    // Function holding all basic matrix operation tests
    bool run2x2MatrixTests();
public:
    SU2Operations();

    bool runSU2Tests();
};

#endif // SU2OPERATIONS_H
