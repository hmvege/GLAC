#ifndef UNITTESTS_H
#define UNITTESTS_H

#include "su3.h"

void SU3BaseTests();
void testOrthogonality(SU3 H, bool verbose);
void testHermicity(SU3 H, bool verbose);
void testNorm(int col, SU3 H);
void testMatrixMultiplication();

#endif // UNITTESTS_H
