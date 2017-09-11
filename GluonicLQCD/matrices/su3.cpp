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
}

SU3::~SU3()
{
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
     */
    complex temp[9];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp[3*i+j] += mat[(3*i+k)]*B[(3*k+j)];
            }
        }
    }
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
        mat[i].re = 0;
        mat[i].im = 0;
    }
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
