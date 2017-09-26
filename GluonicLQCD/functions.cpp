#include "matrices/su3.h"
#include "matrices/su2.h"

void lorentzIndex(int mu, int *lorentzIndices)
{
    /*
     * Fills the lorentz array with correct indices
     */
    for (int i = 0; i < 4; i++)
    {
        lorentzIndices[i] = 0;
    }
    lorentzIndices[mu] = 1;
}

complex SU3Determinant(SU3 U)
{
    complex det;
    det += U[2]*(U[3]*U[7] - U[6]*U[4]);
    det += U[5]*(U[6]*U[1] - U[0]*U[7]);
    det += U[8]*(U[0]*U[4] - U[3]*U[1]);
    return det;
}

SU3 inverse(SU3 U)
{
    /*
     * Gets the inverse of a SU3 matrix, MAKE CLASS BASED?!
     */
    SU3 UInverse; // OPTIMIZE HERE?
    UInverse.copy(U);
    UInverse.conjugate();
    UInverse.transpose();
    return UInverse;
}

bool compareSU3(SU3 A, SU3 B)
{
    for (int i = 0; i < 9; i++) {
        if ((A.mat[i].re() != B.mat[i].re()) || (A.mat[i].im() != B.mat[i].im())) {
            return false;
        }
    }
    return true;
}
