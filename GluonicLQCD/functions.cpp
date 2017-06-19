#include "su3.h"

int index(int i, int j, int k, int l, int N)
{
    /*
     * Function for contigious memory allocation.
     */
    return (N*(N*(N*i + j) + k) + l);
}

SU3 inverse(SU3 U)
{
    SU3 UInverse;
    UInverse.copy(U);
    UInverse.conjugate();
    UInverse.transpose();
    return UInverse;
}
