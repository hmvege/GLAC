#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "matrices/su3.h"
#include "matrices/su2.h"

int index(int i, int j, int k, int l, int N, int N_T);
int indexSpecific(int i, int j, int k, int l, int Ny, int Nz, int Nt);
int stapleIndex(int i, int j, int k, int l, int N, int N_T);
void lorentzIndex(int mu, int *lorentzIndices);
complex SU3Determinant(SU3 U);
SU3 inverse(SU3 U);
SU2 SU2Inverse(SU2 U);
bool compareSU3(SU3 A, SU3 B);

//int neighbourIndex(int proc, int direction, int numprocs);

#endif // FUNCTIONS_H
