/*!
 * \brief File for holding various functions.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "matrices/su3.h"
#include "matrices/su2.h"
#include "complex.h"

inline complex SU2Determinant(const SU2 &H)
{
    // Move into CLASS perhaps best?
    return complex(H.mat[0]*H.mat[6] - H.mat[1]*H.mat[7] - H.mat[2]*H.mat[4] + H.mat[3]*H.mat[5],
            H.mat[0]*H.mat[7] + H.mat[1]*H.mat[6] - H.mat[2]*H.mat[5] - H.mat[3]*H.mat[4]);;
}

inline complex SU3Determinant(const SU3 &H)
{
    /*
     * Function for taking the determinant of
     */
    // TODO: clean up this with a namespace
    return complex( - H.mat[0]*H.mat[10]*H.mat[14] + H.mat[0]*H.mat[11]*H.mat[15] + H.mat[0]*H.mat[16]*H.mat[8] - H.mat[0]*H.mat[17]*H.mat[9] + H.mat[1]*H.mat[10]*H.mat[15] + H.mat[1]*H.mat[11]*H.mat[14]
                    - H.mat[1]*H.mat[16]*H.mat[9] - H.mat[1]*H.mat[17]*H.mat[8] + H.mat[10]*H.mat[12]*H.mat[2] - H.mat[10]*H.mat[13]*H.mat[3] - H.mat[11]*H.mat[12]*H.mat[3] - H.mat[11]*H.mat[13]*H.mat[2]
                    - H.mat[12]*H.mat[4]*H.mat[8] + H.mat[12]*H.mat[5]*H.mat[9] + H.mat[13]*H.mat[4]*H.mat[9] + H.mat[13]*H.mat[5]*H.mat[8] + H.mat[14]*H.mat[4]*H.mat[6] - H.mat[14]*H.mat[5]*H.mat[7]
                    - H.mat[15]*H.mat[4]*H.mat[7] - H.mat[15]*H.mat[5]*H.mat[6] - H.mat[16]*H.mat[2]*H.mat[6] + H.mat[16]*H.mat[3]*H.mat[7] + H.mat[17]*H.mat[2]*H.mat[7] + H.mat[17]*H.mat[3]*H.mat[6],
                    - H.mat[0]*H.mat[10]*H.mat[15] - H.mat[0]*H.mat[11]*H.mat[14] + H.mat[0]*H.mat[16]*H.mat[9] + H.mat[0]*H.mat[17]*H.mat[8] - H.mat[1]*H.mat[10]*H.mat[14] + H.mat[1]*H.mat[11]*H.mat[15]
                    + H.mat[1]*H.mat[16]*H.mat[8] - H.mat[1]*H.mat[17]*H.mat[9] + H.mat[10]*H.mat[12]*H.mat[3] + H.mat[10]*H.mat[13]*H.mat[2] + H.mat[11]*H.mat[12]*H.mat[2] - H.mat[11]*H.mat[13]*H.mat[3]
                    - H.mat[12]*H.mat[4]*H.mat[9] - H.mat[12]*H.mat[5]*H.mat[8] - H.mat[13]*H.mat[4]*H.mat[8] + H.mat[13]*H.mat[5]*H.mat[9] + H.mat[14]*H.mat[4]*H.mat[7] + H.mat[14]*H.mat[5]*H.mat[6]
                    + H.mat[15]*H.mat[4]*H.mat[6] - H.mat[15]*H.mat[5]*H.mat[7] - H.mat[16]*H.mat[2]*H.mat[7] - H.mat[16]*H.mat[3]*H.mat[6] - H.mat[17]*H.mat[2]*H.mat[6] + H.mat[17]*H.mat[3]*H.mat[7]);
}

/*!
 * \brief traceRealMultiplication
 * \param A
 * \param B
 * \return the real trace of the product of the multiplication of A and B.
 */
inline double traceRealMultiplication(const SU3 &A, const SU3 &B)
{
    /*
     * For two regular complex 3x3 matrices trace multiplications taking only real components.
     */
    return (A.mat[0]*B.mat[0] - A.mat[1]*B.mat[1] + A.mat[2]*B.mat[6] - A.mat[3]*B.mat[7] + A.mat[4]*B.mat[12] - A.mat[5]*B.mat[13] +
            A.mat[6]*B.mat[2] - A.mat[7]*B.mat[3] + A.mat[8]*B.mat[8] - A.mat[9]*B.mat[9] + A.mat[10]*B.mat[14] - A.mat[11]*B.mat[15] +
            A.mat[12]*B.mat[4] - A.mat[13]*B.mat[5] + A.mat[14]*B.mat[10] - A.mat[15]*B.mat[11] + A.mat[16]*B.mat[16] - A.mat[17]*B.mat[17]);
}
inline double traceRealMultiplication(SU3 &&A, const SU3 &B)
{
    /*
     * For two regular complex 3x3 matrices trace multiplications taking only real components.
     */
    return (A.mat[0]*B.mat[0] - A.mat[1]*B.mat[1] + A.mat[2]*B.mat[6] - A.mat[3]*B.mat[7] + A.mat[4]*B.mat[12] - A.mat[5]*B.mat[13] +
            A.mat[6]*B.mat[2] - A.mat[7]*B.mat[3] + A.mat[8]*B.mat[8] - A.mat[9]*B.mat[9] + A.mat[10]*B.mat[14] - A.mat[11]*B.mat[15] +
            A.mat[12]*B.mat[4] - A.mat[13]*B.mat[5] + A.mat[14]*B.mat[10] - A.mat[15]*B.mat[11] + A.mat[16]*B.mat[16] - A.mat[17]*B.mat[17]);
}
inline double traceRealMultiplication(const SU3 &A, SU3 &&B)
{
    /*
     * For two regular complex 3x3 matrices trace multiplications taking only real components.
     */
    return (A.mat[0]*B.mat[0] - A.mat[1]*B.mat[1] + A.mat[2]*B.mat[6] - A.mat[3]*B.mat[7] + A.mat[4]*B.mat[12] - A.mat[5]*B.mat[13] +
            A.mat[6]*B.mat[2] - A.mat[7]*B.mat[3] + A.mat[8]*B.mat[8] - A.mat[9]*B.mat[9] + A.mat[10]*B.mat[14] - A.mat[11]*B.mat[15] +
            A.mat[12]*B.mat[4] - A.mat[13]*B.mat[5] + A.mat[14]*B.mat[10] - A.mat[15]*B.mat[11] + A.mat[16]*B.mat[16] - A.mat[17]*B.mat[17]);
}

/*!
 * \brief traceImagMultiplication
 * \param A
 * \param B
 * \return the imaginary trace of the product of the multiplication of A and B.
 */
inline double traceImagMultiplication(const SU3 &A, const SU3 &B)
{
    /*
     * For two regular complex 3x3 matrices trace multiplications taking only imaginary components.
     */
    return (A.mat[0]*B.mat[1] + A.mat[1]*B.mat[0] + A.mat[2]*B.mat[7] + A.mat[3]*B.mat[6] + A.mat[4]*B.mat[13] + A.mat[5]*B.mat[12] +
            A.mat[6]*B.mat[3] + A.mat[7]*B.mat[2] + A.mat[8]*B.mat[9] + A.mat[9]*B.mat[8] + A.mat[10]*B.mat[15] + A.mat[11]*B.mat[14] +
            A.mat[12]*B.mat[5] + A.mat[13]*B.mat[4] + A.mat[14]*B.mat[11] + A.mat[15]*B.mat[10] + A.mat[16]*B.mat[17] + A.mat[17]*B.mat[16]);

}
inline double traceImagMultiplication(SU3 &&A, const SU3 &B)
{
    /*
     * For two regular complex 3x3 matrices trace multiplications taking only imaginary components.
     */
    return (A.mat[0]*B.mat[1] + A.mat[1]*B.mat[0] + A.mat[2]*B.mat[7] + A.mat[3]*B.mat[6] + A.mat[4]*B.mat[13] + A.mat[5]*B.mat[12] +
            A.mat[6]*B.mat[3] + A.mat[7]*B.mat[2] + A.mat[8]*B.mat[9] + A.mat[9]*B.mat[8] + A.mat[10]*B.mat[15] + A.mat[11]*B.mat[14] +
            A.mat[12]*B.mat[5] + A.mat[13]*B.mat[4] + A.mat[14]*B.mat[11] + A.mat[15]*B.mat[10] + A.mat[16]*B.mat[17] + A.mat[17]*B.mat[16]);

}
inline double traceImagMultiplication(const SU3 &A, SU3 &&B)
{
    /*
     * For two regular complex 3x3 matrices trace multiplications taking only imaginary components.
     */
    return (A.mat[0]*B.mat[1] + A.mat[1]*B.mat[0] + A.mat[2]*B.mat[7] + A.mat[3]*B.mat[6] + A.mat[4]*B.mat[13] + A.mat[5]*B.mat[12] +
            A.mat[6]*B.mat[3] + A.mat[7]*B.mat[2] + A.mat[8]*B.mat[9] + A.mat[9]*B.mat[8] + A.mat[10]*B.mat[15] + A.mat[11]*B.mat[14] +
            A.mat[12]*B.mat[5] + A.mat[13]*B.mat[4] + A.mat[14]*B.mat[11] + A.mat[15]*B.mat[10] + A.mat[16]*B.mat[17] + A.mat[17]*B.mat[16]);

}

#endif // FUNCTIONS_H
