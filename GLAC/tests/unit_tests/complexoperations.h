#ifndef COMPLEXOPERATIONS_H
#define COMPLEXOPERATIONS_H

#include "testcore.h"

class ComplexOperations : public TestCore
{
private:
    // Complex class tests
    bool testComplexAddition();
    bool testComplexSubtraction();
    bool testComplexMultiplication();
    bool testComplexDivision();
    bool testComplexConjugation();
    bool testComplexNorm();
    bool testComplexNormSquared();
    bool testComplexSetToMinus();
public:
    ComplexOperations();

    bool runComplexTests();
};

#endif // COMPLEXOPERATIONS_H
