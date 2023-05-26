#ifndef RANDOMMATRIXTESTS_H
#define RANDOMMATRIXTESTS_H

#include "testcore.h"

class RandomMatrixTests : public TestCore
{
private:
    bool testRSTMultiplication();
    bool testRSTInverseMultiplication();
public:
    RandomMatrixTests();

    bool runRandomMatrixTests();
};

#endif // RANDOMMATRIXTESTS_H
