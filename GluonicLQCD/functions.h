#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "su3.h"

int index(int i, int j, int k, int l, int N);
int stapleIndex(int i, int j, int k, int l, int N);
void lorentzIndex(int mu, int *lorentzIndices);
SU3 inverse(SU3 U);

#endif // FUNCTIONS_H
