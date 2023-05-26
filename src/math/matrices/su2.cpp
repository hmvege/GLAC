#include "su2.h"
#include <iostream>
#include <cmath>

/*!
 * \brief SU2::setComplex
 * \param w complex number
 * \param i matrix position(contigious) to set w as.
 */
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
                printf("%12.8f - %12.8fi", mat[4*i+2*j], fabs(mat[4*i+2*j+1]));
            }
            else {
                printf("%12.8f + %12.8fi", mat[4*i+2*j], mat[4*i+2*j+1]);
            }
            printf("      ");
        }
        printf("\n");
    }
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
