#ifndef UNITTESTS_H
#define UNITTESTS_H

// complex tests
void SU3BaseTests();
// Performance tests
void runMatrixPerformanceTest(double eps, double seed, int NTests, bool testMatrix, bool testComplex);


void runBoolTest(int NTests);

#endif // UNITTESTS_H
