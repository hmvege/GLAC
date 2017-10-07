#include "su3.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "complex.h"

using std::cout;
using std::endl;

SU3::SU3()
{
    // Remove this! Make user manually fill by zeros! Can save quite some time (y)
//    for (int i = 0; i < 18; i++) {
//        mat[i] = 0;
//    }
}

SU3::~SU3()
{
}

SU3 &SU3::operator=(const SU3 &B)
{
    for (int i = 0; i < 18; i++) {
        mat[i] = B.mat[i];
    }
    return *this;
}

SU3 &SU3::operator+=(SU3 B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

SU3 &SU3::operator-=(SU3 B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

SU3 &SU3::operator*=(SU3 B)
{
    /*
     * ab = (a + bi)(c + id) = ac + iad + ibc - bd = ac - bd + i(ad + bc);
     * 0 1 2   11 12 13   00 01 02
     * 3 4 5 = 21 22 23 = 10 11 12
     * 6 7 8   31 32 33   20 21 22
     */
    double temp[18];

    temp[0] = mat[0]*B[0] - mat[1]*B[1] + mat[2]*B[6] - mat[3]*B[7] + mat[4]*B[12] - mat[5]*B[13];
    temp[1] = mat[0]*B[1] + mat[1]*B[0] + mat[2]*B[7] + mat[3]*B[6] + mat[4]*B[13] + mat[5]*B[12];
    temp[2] = mat[0]*B[2] - mat[1]*B[3] + mat[2]*B[8] - mat[3]*B[9] + mat[4]*B[14] - mat[5]*B[15];
    temp[3] = mat[0]*B[3] + mat[1]*B[2] + mat[2]*B[9] + mat[3]*B[8] + mat[4]*B[15] + mat[5]*B[14];
    temp[4] = mat[0]*B[4] - mat[1]*B[5] + mat[2]*B[10] - mat[3]*B[11] + mat[4]*B[16] - mat[5]*B[17];
    temp[5] = mat[0]*B[5] + mat[1]*B[4] + mat[2]*B[11] + mat[3]*B[10] + mat[4]*B[17] + mat[5]*B[16];
    temp[6] = mat[6]*B[0] - mat[7]*B[1] + mat[8]*B[6] - mat[9]*B[7] + mat[10]*B[12] - mat[11]*B[13];
    temp[7] = mat[6]*B[1] + mat[7]*B[0] + mat[8]*B[7] + mat[9]*B[6] + mat[10]*B[13] + mat[11]*B[12];
    temp[8] = mat[6]*B[2] - mat[7]*B[3] + mat[8]*B[8] - mat[9]*B[9] + mat[10]*B[14] - mat[11]*B[15];
    temp[9] = mat[6]*B[3] + mat[7]*B[2] + mat[8]*B[9] + mat[9]*B[8] + mat[10]*B[15] + mat[11]*B[14];
    temp[10] = mat[6]*B[4] - mat[7]*B[5] + mat[8]*B[10] - mat[9]*B[11] + mat[10]*B[16] - mat[11]*B[17];
    temp[11] = mat[6]*B[5] + mat[7]*B[4] + mat[8]*B[11] + mat[9]*B[10] + mat[10]*B[17] + mat[11]*B[16];
    temp[12] = mat[12]*B[0] - mat[13]*B[1] + mat[14]*B[6] - mat[15]*B[7] + mat[16]*B[12] - mat[17]*B[13];
    temp[13] = mat[12]*B[1] + mat[13]*B[0] + mat[14]*B[7] + mat[15]*B[6] + mat[16]*B[13] + mat[17]*B[12];
    temp[14] = mat[12]*B[2] - mat[13]*B[3] + mat[14]*B[8] - mat[15]*B[9] + mat[16]*B[14] - mat[17]*B[15];
    temp[15] = mat[12]*B[3] + mat[13]*B[2] + mat[14]*B[9] + mat[15]*B[8] + mat[16]*B[15] + mat[17]*B[14];
    temp[16] = mat[12]*B[4] - mat[13]*B[5] + mat[14]*B[10] - mat[15]*B[11] + mat[16]*B[16] - mat[17]*B[17];
    temp[17] = mat[12]*B[5] + mat[13]*B[4] + mat[14]*B[11] + mat[15]*B[10] + mat[16]*B[17] + mat[17]*B[16];

    for (int i = 0; i < 18; i++) // OPTIMIZATION: can write the last element directly to the mat[18] instead of through temp
    {
        mat[i] = temp[i];
    }
    return *this;
}

SU3 &SU3::operator/=(double a)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] /= a;
    }
    return *this;
}

SU3 SU3::inv()
{   /*
     * Takes the inverse of the matrix(which is transpose and conjugate).
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
    SU3 R;
    // Upper triangular
    R.mat[6]  =  mat[2];
    R.mat[7]  = -mat[3];
    R.mat[12] =  mat[4];
    R.mat[13] = -mat[5];
    R.mat[14] =  mat[10];
    R.mat[15] = -mat[11];
    // Upper triangular
    R.mat[2]  =  mat[6];
    R.mat[3]  = -mat[7];
    R.mat[4]  =  mat[12];
    R.mat[5]  = -mat[13];
    R.mat[10] =  mat[14];
    R.mat[11] = -mat[15];
    // Diagonals
    R.mat[0]  =  mat[0];
    R.mat[1]  = -mat[1];
    R.mat[8]  =  mat[8];
    R.mat[9]  = -mat[9];
    R.mat[16] =  mat[16];
    R.mat[17] = -mat[17];
    return R;
}

void SU3::zeros()
{
    for (int i = 0; i < 18; i++) {
        mat[i] = 0;
    }
}

void SU3::identity()
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] = 0;
    }
    mat[0] = 1;
    mat[8] = 1;
    mat[16] = 1;
}

SU3 SU3::transpose()
{
    double temp[18];
    for (int i = 0; i < 18; i++) { temp[i] = mat[i]; }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                mat[i*6+2*j+k] = temp[j*6+2*i+k];
                mat[j*6+2*i+k] = temp[i*6+2*j+k];
            }
        }
    }
    return *this;
}

SU3 SU3::conjugate()
{
    for (int i = 0; i < 9; i++) {
        mat[2*i+1] *= -1;
    }
    return *this;
}

void SU3::copy(SU3 B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] = B.mat[i];
    }
}

double SU3::norm(int i)
{
    /*
     * Returns the norm of the complex number.
     */
    return sqrt(mat[i]*mat[i] + mat[i+1]*mat[i+1]);
}

double SU3::normSquared(int i)
{
    /*
     * Returns the norm squared of the complex number.
     */
    return mat[i]*mat[i] + mat[i+1]*mat[i+1];
}

complex SU3::c(int i)
{
    /*
     * Returns the complex conjugate of the object instance.
     */
//    mat[i+1] *= -1;
//    return *this;
    return complex(mat[2*i],mat[2*i+1]);
}

//complex SU3::getComplex(int i)
//{
//    return complex(mat[2*i],mat[2*i+1]);
//}

void SU3::setComplex(complex w, int i)
{
    mat[i] = w.z[0];
    mat[i+1] = w.z[1];
}

void SU3::print()
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout <<std::setw(12) << mat[6*i + 2*j];
            if (mat[6*i + 2*j + 1] < 0)
            {
                cout << " - ";
            }
            else
            {
                cout << " + ";
            }
            cout << std::setw(12) << fabs(mat[6*i + 2*j + 1]) << "i,  ";
        }
        cout << endl;
    }
}
