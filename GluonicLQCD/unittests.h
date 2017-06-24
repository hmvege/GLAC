#ifndef UNITTESTS_H
#define UNITTESTS_H

#include "matrices/su3.h"

void SU3BaseTests();
bool testOrthogonality(SU3 H, bool verbose);
bool testHermicity(SU3 H, bool verbose);
bool testNorm(int col, SU3 H);
void testMatrixMultiplication();
bool SU2UnitTest(complex * r, complex * s, complex * t);

#endif // UNITTESTS_H
