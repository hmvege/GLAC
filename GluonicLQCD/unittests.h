#ifndef UNITTESTS_H
#define UNITTESTS_H

#include "su3.h"

void SU3BaseTests();
bool testOrthogonality(SU3 H, bool verbose);
bool testHermicity(SU3 H, bool verbose);
bool testNorm(int col, SU3 H);
void testMatrixMultiplication();

#endif // UNITTESTS_H
