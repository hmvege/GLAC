#include "su3.h"
#include "complex.h"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

/*
 * For storing the SU3 matrix
 */

SU3::SU3()
{
//    mat = new complex[9];
}

SU3::~SU3()
{
//    delete [] mat;
}

void SU3::transpose()
{
    complex temp[9];
    for (int i = 0; i < 9; i++) { temp[i] = mat[i]; }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < i; j++)
        {
            mat[i*3+j] = temp[j*3+i];
            mat[j*3+i] = temp[i*3+j];
        }
    }
}

void SU3::conjugate()
{
    for (int i = 0; i < 9; i++)
    {
        mat[i].conjugate();
    }
}

void SU3::copy(SU3 B)
{
    for (int i = 0; i < 9; i++)
    {
        mat[i] = B.mat[i];
    }
}

SU3 SU3::inv()
{
    SU3 R;
    // Upper triangular
    R.mat[3] = mat[1].c();
    R.mat[6] = mat[2].c();
    R.mat[7] = mat[5].c();
    // Upper triangular
    R.mat[1] = mat[3].c();
    R.mat[2] = mat[6].c();
    R.mat[5] = mat[7].c();
    // Diagonals
    R.mat[0] = mat[0].c();
    R.mat[4] = mat[4].c();
    R.mat[8] = mat[8].c();
    return R;
}

SU3 &SU3::operator=(const SU3 &B)
{
    for (int i = 0; i < 9; i++) {
        mat[i] = B.mat[i];
    }
    return *this;
}

SU3 &SU3::operator+=(SU3 B)
{
    for (int i = 0; i < 9; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

SU3 &SU3::operator*=(SU3 B)
{
    /*
     * a*b = (a + bi)(c + id) = a*c + iad + ibc - bd;
     * 0 1 2   11 12 13
     * 3 4 5 = 21 22 23
     * 6 7 8   31 32 33
     */
    complex temp[9];
    temp[0] = mat[0]*B.mat[0] + mat[1]*B.mat[3] + mat[2]*B.mat[6];
    temp[1] = mat[0]*B.mat[1] + mat[1]*B.mat[4] + mat[2]*B.mat[7];
    temp[2] = mat[0]*B.mat[2] + mat[1]*B.mat[5] + mat[2]*B.mat[8];
    temp[3] = mat[3]*B.mat[0] + mat[4]*B.mat[3] + mat[5]*B.mat[6];
    temp[4] = mat[3]*B.mat[1] + mat[4]*B.mat[4] + mat[5]*B.mat[7];
    temp[5] = mat[3]*B.mat[2] + mat[4]*B.mat[5] + mat[5]*B.mat[8];
    temp[6] = mat[6]*B.mat[0] + mat[7]*B.mat[3] + mat[8]*B.mat[6];
    temp[7] = mat[6]*B.mat[1] + mat[7]*B.mat[4] + mat[8]*B.mat[7];
    temp[8] = mat[6]*B.mat[2] + mat[7]*B.mat[5] + mat[8]*B.mat[8];
//    for (int i = 0; i < 3; i++) // This is the fastest method for some reason
//    {
//        for (int j = 0; j < 3; j++)
//        {
//            for (int k = 0; k < 3; k++)
//            {
//                temp[3*i+j] += mat[(3*i+k)]*B[(3*k+j)];
//            }
//        }
//    }
    for (int i = 0; i < 9; i++)
    {
        mat[i] = temp[i];
    }
    return *this;
}

SU3 &SU3::operator-=(SU3 B)
{
    for (int i = 0; i < 9; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

void SU3::zeros()
{
    for (int i = 0; i < 9; i++)
    {
        mat[i].setRe(0);
        mat[i].setIm(0);
    }
}

void SU3::identity()
{
    for (int i = 0; i < 9; i++)
    {
        mat[i].setRe(0);
        mat[i].setIm(0);
    }
    mat[0].setRe(1);
    mat[4].setRe(1);
    mat[8].setRe(1);
}



void SU3::print()
{
    /*
     * Temporary class for printing matrix. Might remove in the future to get better performance
     */
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << std::setw(15) << mat[3*i+j];
        }
        cout << endl;
    }
}
