#include "su2.h"

SU2::SU2()
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = 0;
    }
}

SU2::~SU2()
{
    /*
     * Destructor operator
     */
}

SU2 &SU2::operator=(const SU2 &B)
{
    /*
     * Copy assignement operator.
     */
    for (int i = 0; i < 8; i++) {
        mat[i] = B.mat[i];
    }
    return *this;
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

    temp[0] = mat[0]*B[0] - mat[1]*B[1] + mat[2]*B[4] - mat[3]*B[5];
    temp[1] = mat[0]*B[1] + mat[1]*B[0] + mat[2]*B[5] + mat[3]*B[4];
    temp[2] = mat[0]*B[2] - mat[1]*B[3] + mat[2]*B[6] - mat[3]*B[7];
    temp[3] = mat[0]*B[3] + mat[1]*B[2] + mat[2]*B[7] + mat[3]*B[6];
    temp[4] = mat[4]*B[0] - mat[5]*B[1] + mat[6]*B[4] - mat[7]*B[5];
    temp[5] = mat[4]*B[1] + mat[5]*B[0] + mat[6]*B[5] + mat[7]*B[4];
    temp[6] = mat[4]*B[2] - mat[5]*B[3] + mat[6]*B[6] - mat[7]*B[7];
    temp[7] = mat[4]*B[3] + mat[5]*B[2] + mat[6]*B[7] + mat[7]*B[6];

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

void SU2::setComplex(complex w, int i)
{
    mat[i] = w.z[0];
    mat[i+1] = w.z[1];
}

double SU2::normSquared(int i)
{
    /*
     * Returns the norm squared of the complex number.
     */
    return mat[i]*mat[i] + mat[i+1]*mat[i+1];
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

SU2 SU2::inv()
{
    SU2 R;
    R[0] = mat[0];
    R[1] = -mat[1];

    R[2] = mat[4];
    R[3] = -mat[5];

    R[4] = mat[2];
    R[5] = -mat[3];

    R[6] = mat[6];
    R[7] = -mat[7];

    return R;
}

SU2 SU2::transpose()
{
    double temp[8];
    for (int i = 0; i < 8; i++) { temp[i] = mat[i]; }

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                mat[i*4+2*j+k] = temp[j*4+2*i+k];
                mat[j*4+2*i+k] = temp[i*4+2*j+k];
            }
        }
    }
    return *this;
}

SU2 SU2::conjugate()
{
    for (int i = 0; i < 4; i++) {
        mat[2*i+1] = -mat[2*i+1];
    }
    return *this;
}
