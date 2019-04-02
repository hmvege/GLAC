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

    /*
     * identity, su3, complex, double
     * zeros, su3, complex, double
     * sumSpatial, doubles
     * realTraceMultiplication, SU3
     * imagTraceMultiplication, SU3
     * transpose
     * conjugate
     * makeAntiHermitian
     * makeHermitian
     *
     */

    bool testMakeHermitian();
    bool testMakeAntiHermitian();
    bool testLatticeSumSpatial();
    bool testLatticeRealTraceMultiplication();
    bool testLatticeImagTraceMultiplication();
    bool testTranpose();
    bool testConjugate();
    bool testIdentity();
    bool testZeros();

    // Tests including parallel communication
    bool fullLatticeTests();
    bool testLatticeShift();
    bool testFieldGaugeInvariance();
public:
    LatticeOperations();

    bool runLatticeTests();
};

#endif // LATTICEOPERATIONS_H
