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
    /*
     * Gets the inverse of a SU3 matrix, MAKE CLASS BASED?!
     */
    SU3 UInverse;
    UInverse.copy(U);
    UInverse.conjugate();
    UInverse.transpose();
    return UInverse;
}

int stapleIndex(int i, int j, int k, int l, int N)
{
    /*
     * Unit vector for lorentz indexes
     */
    return index((i+N) % N, (j+N) % N, (k+N) % N, (l+N) % N, N);;
}

void lorentzIndex(int mu, int *lorentzIndices)
{
    /*
     * Fills the lorentz array with correct indices
     */
    for (int i = 0; i < 4; i++)
    {
        if (mu==i)
        {
            lorentzIndices[i] = 1;
        }
        else
        {
            lorentzIndices[i] = 0;
        }
    }
}
