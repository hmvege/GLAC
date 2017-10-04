#include "_su3.h"

_su3::_su3()
{
    for (int i = 0; i < 18; i++) {
        mat[i] = 0;
    }
}

_su3::~_su3()
{
}

_su3 &_su3::operator=(const _su3 &B)
{
    for (int i = 0; i < 18; i++) {
        mat[i] = B.mat[i];
    }
    return *this;
}

_su3 &_su3::operator+=(_su3 B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

_su3 &_su3::operator-=(_su3 B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

_su3 &_su3::operator*=(_su3 B)
{
    /*
     * a*b = (a + bi)(c + id) = a*c + iad + ibc - bd;
     * 0 1 2   11 12 13   00 01 02
     * 3 4 5 = 21 22 23 = 10 11 12
     * 6 7 8   31 32 33   20 21 22
     */
    double temp[18];

    temp[0] = mat[0]*B[0] - mat[1]*B[1] + mat[0]*B[2] - mat[1]*B[3] + mat[0]*B[4] - mat[1]*B[5];
    temp[1] = mat[0]*B[1] + mat[1]*B[0] + mat[0]*B[3] + mat[1]*B[2] + mat[0]*B[5] + mat[1]*B[4];
    temp[2] = mat[0]*B[0] - mat[1]*B[1] + mat[0]*B[2] - mat[1]*B[3] + mat[0]*B[4] - mat[1]*B[5];
    temp[3] = mat[0]*B[1] + mat[1]*B[0] + mat[0]*B[3] + mat[1]*B[2] + mat[0]*B[5] + mat[1]*B[4];
    temp[4] = mat[0]*B[0] - mat[1]*B[1] + mat[0]*B[2] - mat[1]*B[3] + mat[0]*B[4] - mat[1]*B[5];
    temp[5] = mat[0]*B[1] + mat[1]*B[0] + mat[0]*B[3] + mat[1]*B[2] + mat[0]*B[5] + mat[1]*B[4];
    temp[6] = mat[2]*B[0] - mat[3]*B[1] + mat[2]*B[2] - mat[3]*B[3] + mat[2]*B[4] - mat[3]*B[5];
    temp[7] = mat[2]*B[1] + mat[3]*B[0] + mat[2]*B[3] + mat[3]*B[2] + mat[2]*B[5] + mat[3]*B[4];
    temp[8] = mat[2]*B[0] - mat[3]*B[1] + mat[2]*B[2] - mat[3]*B[3] + mat[2]*B[4] - mat[3]*B[5];
    temp[9] = mat[2]*B[1] + mat[3]*B[0] + mat[2]*B[3] + mat[3]*B[2] + mat[2]*B[5] + mat[3]*B[4];
    temp[10] = mat[2]*B[0] - mat[3]*B[1] + mat[2]*B[2] - mat[3]*B[3] + mat[2]*B[4] - mat[3]*B[5];
    temp[11] = mat[2]*B[1] + mat[3]*B[0] + mat[2]*B[3] + mat[3]*B[2] + mat[2]*B[5] + mat[3]*B[4];
    temp[12] = mat[4]*B[0] - mat[5]*B[1] + mat[4]*B[2] - mat[5]*B[3] + mat[4]*B[4] - mat[5]*B[5];
    temp[13] = mat[4]*B[1] + mat[5]*B[0] + mat[4]*B[3] + mat[5]*B[2] + mat[4]*B[5] + mat[5]*B[4];
    temp[14] = mat[4]*B[0] - mat[5]*B[1] + mat[4]*B[2] - mat[5]*B[3] + mat[4]*B[4] - mat[5]*B[5];
    temp[15] = mat[4]*B[1] + mat[5]*B[0] + mat[4]*B[3] + mat[5]*B[2] + mat[4]*B[5] + mat[5]*B[4];
    temp[16] = mat[4]*B[0] - mat[5]*B[1] + mat[4]*B[2] - mat[5]*B[3] + mat[4]*B[4] - mat[5]*B[5];
    temp[17] = mat[4]*B[1] + mat[5]*B[0] + mat[4]*B[3] + mat[5]*B[2] + mat[4]*B[5] + mat[5]*B[4];

    for (int i = 0; i < 18; i++)
    {
        mat[i] = temp[i];
    }
    return *this;
}

_su3 _su3::inv()
{
    _su3 R;
    // Upper triangular
    R.mat[6] = mat[2];
    R.mat[12] = mat[4];
    R.mat[14] = mat[10];
    R.mat[6] = mat[3];
    R.mat[12] = mat[5];
    R.mat[14] = mat[11];
    // Upper triangular
    R.mat[2] = mat[6];
    R.mat[3] = mat[7];
    R.mat[4] = mat[12];
    R.mat[5] = mat[13];
    R.mat[10] = mat[14];
    R.mat[11] = mat[15];
    // Diagonals
    R.mat[0]  =  mat[0];
    R.mat[1]  = -mat[1];
    R.mat[8]  =  mat[8];
    R.mat[9]  = -mat[9];
    R.mat[16] =  mat[16];
    R.mat[17] = -mat[17];
    return R;
}

void _su3::zeros()
{
    for (int i = 0; i < 18; i++) {
        mat[i] = 0;
    }
}

void _su3::identity()
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] = 0;
    }
    mat[0] = 1;
    mat[8] = 1;
    mat[16] = 1;
}

void _su3::transpose()
{
    double temp[18];
    for (int i = 0; i < 18; i++) { temp[i] = mat[i]; }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                mat[i*3+j+k] = temp[j*3+i+k];
                mat[j*3+i+k] = temp[i*3+j+k];
            }
        }
    }
}

void _su3::conjugate()
{
    for (int i = 0; i < 9; i++) {
        mat[2*i+1] = -mat[2*i+1];
    }
}

void _su3::copy(_su3 B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] = B.mat[i];
    }
}
