#include "su3.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

///////////////////////////////
//// SU3 specific functions ///
///////////////////////////////

SU3 SU3::inv()
{
    /*
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
    return complex(mat[2*i],mat[2*i+1]);
}

void SU3::setComplex(complex w, int i)
{
    mat[i] = w.z[0];
    mat[i+1] = w.z[1];
}

complex SU3::trace()
{
    return complex(mat[0] + mat[8] + mat[16], mat[1] + mat[9] + mat[17]);
}

SU3 SU3::makeHermitian()
{
    /*
     * An anti-hermitian matrix is made hermitian by multiplying by (-i)
     */
    double temp;
    for (int i = 0; i < 9; i++) {
        temp = mat[2*i];
        mat[2*i] = mat[2*i+1];
        mat[2*i+1] = -temp;
    }
    return *this;
}

SU3 SU3::makeAntiHermitian()
{
    /*
     * Multiplies by (i). Ensure this is correct in unit tests!
     */
    double temp;
    for (int i = 0; i < 9; i++) {
        temp = mat[2*i];
        mat[2*i] = -mat[2*i+1];
        mat[2*i+1] = temp;
    }
    return *this;
}

SU3 SU3::getIm()
{
    for (int i = 0; i < 18; i+=2)
    {
        mat[i] = 0;
    }
    return *this;
}

SU3 SU3::getRe()
{
    for (int i = 1; i < 18; i+=2)
    {
        mat[i] = 0;
    }
    return *this;
}

void SU3::print()
{
    for (int i = 0; i < 3; i++) // Machine friendly way
    {
        if (i == 0) {
            printf("[");
        } else {
            printf(" ");
        }
        printf("[");
        for (int j = 0; j < 3; j++)
        {
            printf("%12.8f", mat[6*i + 2*j]);
            if (mat[6*i + 2*j + 1] < 0)
            {
                printf(" - ");
            }
            else
            {
                printf(" + ");
            }
            printf("%12.8fi",fabs(mat[6*i + 2*j + 1]));
            if (j == 2){
                printf("]");
                if (i != 2) printf(",");
            } else {
                printf(", ");
            }

        }
        if (i != 2) printf("\n");
    }
    printf("]\n");
}
