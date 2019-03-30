#ifndef LATTICEOPERATIONS_H
#define LATTICEOPERATIONS_H

#include "testcore.h"

class LatticeOperations : public TestCore
{
private:
    // Lattice class tests
    bool testLatticeAddition();
    bool testLatticeSubtraction();
    bool testLatticeMultiplication();
    bool testLatticeDivision();
    bool testLatticeTrace();
    bool testLatticeRealTrace();
    bool testLatticeImagTrace();
    bool testLatticeSubtractReal();
    bool testLatticeSubtractImag();
    bool testLatticeSum();
    bool testLatticeSumRealTrace();
    bool testLatticeSumRealTraceMultiplication();
    bool testLatticeInverse();

    // Tests including parallel communication
    bool fullLatticeTests();
    bool testLatticeShift();
    bool testFieldGaugeInvariance();
public:
    LatticeOperations();

    bool runLatticeTests();
};

#endif // LATTICEOPERATIONS_H
