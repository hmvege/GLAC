#include "su2.h"

using std::cout;
using std::endl;

SU2::SU2()
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = 0;
    }
}

SU2::~SU2()
{
}

SU2 &SU2::operator+=(SU2 B)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

SU2 &SU2::operator-=(SU2 B)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

SU2 &SU2::operator*=(double b)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] *= b;
    }
    return *this;
}

SU2 &SU2::operator*=(SU2 B)
{
    double temp[8];

    temp[0] = mat[0]*B[0] - mat[1]*B[1] + mat[0]*B[2] - mat[1]*B[3];
    temp[1] = mat[0]*B[1] + mat[1]*B[0] + mat[0]*B[3] + mat[1]*B[2];
    temp[2] = mat[0]*B[0] - mat[1]*B[1] + mat[0]*B[2] - mat[1]*B[3];
    temp[3] = mat[0]*B[1] + mat[1]*B[0] + mat[0]*B[3] + mat[1]*B[2];
    temp[4] = mat[2]*B[0] - mat[3]*B[1] + mat[2]*B[2] - mat[3]*B[3];
    temp[5] = mat[2]*B[1] + mat[3]*B[0] + mat[2]*B[3] + mat[3]*B[2];
    temp[6] = mat[2]*B[0] - mat[3]*B[1] + mat[2]*B[2] - mat[3]*B[3];
    temp[7] = mat[2]*B[1] + mat[3]*B[0] + mat[2]*B[3] + mat[3]*B[2];

    for (int i = 0; i < 8; i++)
    {
        mat[i] = temp[i];
    }
    return *this;
}

void SU2::zeros()
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = 0;
    }
}

void SU2::identity()
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = 0;
    }
    mat[0] = 1;
    mat[6] = 1;
}

void SU2::copy(SU2 B)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = B.mat[i];
    }
}

void SU2::print()
{
    /*
     * Temporary class for printing matrix. Might remove in the future to get better performance
     */
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if (mat[4*i+2*j+1] < 0) {
                cout << std::setw(15) << mat[4*i+2*j] << " - " << fabs(mat[4*i+2*j+1]) << "i";
            }
            else {
                cout << std::setw(15) << mat[4*i+2*j] << " + " << mat[4*i+2*j+1] << "i";
            }
            cout << "  ";
        }
        cout << endl;
    }
}
