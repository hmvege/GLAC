#include "su3.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;

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

complex SU3::trace()
{
    return complex(mat[0] + mat[8] + mat[16], mat[1] + mat[9] + mat[17]);
}

void SU3::setComplex(complex w, int i)
{
    mat[i] = w.z[0];
    mat[i+1] = w.z[1];
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

// Overloading stream operator
std::ostream& operator<<(std::ostream& os, const SU3& mat)
{
    const auto prec = std::cout.precision();
    os << std::setprecision(2);
    os << std::scientific;

    for (int i = 0; i < 3; i++)
    {
        if (i == 0) {
            os << "[";
        } else {
            os << " ";
        }
        os << "[";
        for (int j = 0; j < 3; j++)
        {
            os << std::setw(8) << mat.mat[6*i + 2*j];
            if (mat.mat[6*i + 2*j + 1] < 0)
            {
                os << " - ";
            }
            else
            {
                os << " + ";
            }
            os << std::setw(8) << std::fabs(mat.mat[6*i + 2*j + 1]) << "i";
            if (j == 2)
            {
                {
                    os << "]";
                }
                if (i != 2)
                {
                    os << ",";
                }
            }
            else
            {
                os << ", ";
            }

        }
        if (i != 2)
        {
            os << "\n";
        }
    }
    os << "]";
    os << std::setprecision(prec);

    return os;
}
