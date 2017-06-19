#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "su3.h"

int index(int i, int j, int k, int l, int N);
int indexAction(int i, int j, int k, int l, int mu, int N);
SU3 inverse(SU3 U);

#endif // FUNCTIONS_H
