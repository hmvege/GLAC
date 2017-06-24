#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "matrices/su3.h"
#include "matrices/su2.h"

int index(int i, int j, int k, int l, int N);
int stapleIndex(int i, int j, int k, int l, int N, int N_T);
void lorentzIndex(int mu, int *lorentzIndices);
SU3 inverse(SU3 U);
SU2 SU2Inverse(SU2 U);

#endif // FUNCTIONS_H
