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
    for (int i = 0; i < 9; i++)
    {
        mat[i].re = 0;
        mat[i].im = 0; // TEMP set to i+1
    }
}

SU3::~SU3()
{

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
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            fortsetther!
            cout << "TODO: FIX MATRIX MULTIPLICATION" << endl;
        }
    }
    return *this;
}

SU3 &SU3::operator-=(SU3 B)
{
    for (int i = 0; i < 9; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
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
            cout << std::setw(10) << mat[3*i+j];
        }
        cout << endl;
    }
}
