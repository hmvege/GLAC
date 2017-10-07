#include "matrices/su3.h"
#include "matrices/su2.h"
#include "complex.h"

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
    // Redo! Move into class perhaps best?
    complex det;
    det.setRe(U[2]*(U[6]*U[14] - U[7]*U[15] - U[12]*U[8] + U[13]*U[9]) - U[3]*(U[6]*U[15] + U[7]*U[14] - U[12]*U[9] - U[13]*U[8])
            + U[10]*(U[12]*U[2] - U[13]*U[3] - (U[0]*U[14] - U[1]*U[15])) - U[11]*(U[12]*U[3] + U[13]*U[2] - (U[0]*U[15] + U[1]*U[14]))
            + U[16]*(U[0]*U[8] - U[1]*U[9] - (U[6]*U[2] - U[7]*U[3])) - U[17]*(U[0]*U[9] + U[1]*U[8] - (U[6]*U[3] + U[7]*U[2])));
    det.setIm(U[2]*(U[6]*U[15] + U[7]*U[14] - U[12]*U[9] - U[13]*U[8]) + U[3]*(U[6]*U[14] - U[7]*U[15] - U[12]*U[8] + U[13]*U[9])
            + U[10]*(U[12]*U[3] + U[13]*U[2] - (U[0]*U[15] + U[1]*U[14])) - U[11]*(U[12]*U[2] - U[13]*U[3] - (U[0]*U[14] - U[1]*U[15]))
            + U[16]*(U[0]*U[9] + U[1]*U[8] - (U[6]*U[3] + U[7]*U[2])) - U[17]*(U[0]*U[8] - U[1]*U[9] - (U[6]*U[2] - U[7]*U[3])));
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

//bool compareSU3(SU3 A, SU3 B)
//{
//    // Redo this as well... with epsilon?
//    for (int i = 0; i < 18; i++) {
//        if (A.mat[i] != B.mat[i]) {
//            return false;
//        }
//    }
//    return true;
//}

double traceRealMultiplication(SU3 A, SU3 B)
{
    // Redo this in python
    return (A[0]*B[0] - A[1]*B[1] + A[2]*B[6] - A[3]*B[7] + A[4]*B[12] - A[5]*B[13] +
            A[6]*B[2] - A[7]*B[3] + A[8]*B[8] - A[9]*B[9] + A[10]*B[14] - A[11]*B[15] +
            A[12]*B[4] - A[13]*B[5] + A[14]*B[10] - A[15]*B[11] + A[16]*B[16] - A[17]*B[17]);

//    return ((A.mat[0]*B.mat[0] + A.mat[1]*B.mat[3] + A.mat[2]*B.mat[6]).re() + // COMPLEX MULT HERE
//            (A.mat[3]*B.mat[1] + A.mat[4]*B.mat[4] + A.mat[5]*B.mat[7]).re() +
//            (A.mat[6]*B.mat[2] + A.mat[7]*B.mat[5] + A.mat[8]*B.mat[8]).re());
}

complex complexMultiply(SU3 A, SU3 B, int i, int j)
{
    /*
     * For multiplying two complex numbers of the new SU3 matrix classes
     * a*b = (a + bi)(c + id) = ac + iad + ibc - bd = ac - bd + i(ad + bc);
     */
    return complex(A[i]*B[j] - A[i+1]*B[j+1],
                   A[i]*B[j+1] + A[i+1]*B[j]);
}
