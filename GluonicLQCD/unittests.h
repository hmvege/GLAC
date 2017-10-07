#ifndef UNITTESTS_H
#define UNITTESTS_H

#include "matrices/su3.h"

// complex tests
void SU3BaseTests();
// Performance tests
void runMatrixPerformanceTest(double eps, double seed, int NTests, bool testMatrix, bool testComplex);


#endif // UNITTESTS_H
