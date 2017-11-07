#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "matrices/su3.h"
#include "matrices/su2.h"
#include "complex.h"

complex SU2Determinant(SU2 H);
complex SU3Determinant(SU3 U);
double traceRealMultiplication(SU3 A, SU3 B);
double traceImagMultiplication(SU3 A, SU3 B);
double traceSparseRealMultiplication(SU3 A, SU3 B);
double traceSparseImagMultiplication(SU3 A, SU3 B);
complex complexMultiply(SU3 A, SU3 B, int i, int j);

#endif // FUNCTIONS_H
