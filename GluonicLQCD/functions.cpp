#include "matrices/su3.h"
#include "matrices/su2.h"
#include "complex.h"

complex SU3Determinant(SU3 H)
{
    // Redo! Move into class perhaps best?
    return complex( - H[0]*H[10]*H[14] + H[0]*H[11]*H[15] + H[0]*H[16]*H[8] - H[0]*H[17]*H[9] + H[1]*H[10]*H[15] + H[1]*H[11]*H[14]
                    - H[1]*H[16]*H[9] - H[1]*H[17]*H[8] + H[10]*H[12]*H[2] - H[10]*H[13]*H[3] - H[11]*H[12]*H[3] - H[11]*H[13]*H[2]
                    - H[12]*H[4]*H[8] + H[12]*H[5]*H[9] + H[13]*H[4]*H[9] + H[13]*H[5]*H[8] + H[14]*H[4]*H[6] - H[14]*H[5]*H[7]
                    - H[15]*H[4]*H[7] - H[15]*H[5]*H[6] - H[16]*H[2]*H[6] + H[16]*H[3]*H[7] + H[17]*H[2]*H[7] + H[17]*H[3]*H[6],
                    - H[0]*H[10]*H[15] - H[0]*H[11]*H[14] + H[0]*H[16]*H[9] + H[0]*H[17]*H[8] - H[1]*H[10]*H[14] + H[1]*H[11]*H[15]
                    + H[1]*H[16]*H[8] - H[1]*H[17]*H[9] + H[10]*H[12]*H[3] + H[10]*H[13]*H[2] + H[11]*H[12]*H[2] - H[11]*H[13]*H[3]
                    - H[12]*H[4]*H[9] - H[12]*H[5]*H[8] - H[13]*H[4]*H[8] + H[13]*H[5]*H[9] + H[14]*H[4]*H[7] + H[14]*H[5]*H[6]
                    + H[15]*H[4]*H[6] - H[15]*H[5]*H[7] - H[16]*H[2]*H[7] - H[16]*H[3]*H[6] - H[17]*H[2]*H[6] + H[17]*H[3]*H[7]);
}

double traceRealMultiplication(SU3 A, SU3 B)
{
    // Redo this in python
    return (A[0]*B[0] - A[1]*B[1] + A[2]*B[6] - A[3]*B[7] + A[4]*B[12] - A[5]*B[13] +
            A[6]*B[2] - A[7]*B[3] + A[8]*B[8] - A[9]*B[9] + A[10]*B[14] - A[11]*B[15] +
            A[12]*B[4] - A[13]*B[5] + A[14]*B[10] - A[15]*B[11] + A[16]*B[16] - A[17]*B[17]);
}

complex complexMultiply(SU3 A, SU3 B, int i, int j)
{
    /*
     * For multiplying two complex numbers of the new SU3 matrix classes
     * a*b = (a + bi)(c + id) = ac + iad + ibc - bd = ac - bd + i(ad + bc);
     */
    std::cout << "This wierd complexMultiply in functions.cpp is used!" <<  std::endl;
    return complex(A[i]*B[j] - A[i+1]*B[j+1],
                   A[i]*B[j+1] + A[i+1]*B[j]);
}
