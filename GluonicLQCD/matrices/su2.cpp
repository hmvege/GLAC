#include "su2.h"


#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;


SU2::SU2()
{
}

SU2::~SU2()
{
}

void SU2::transpose()
{
    complex temp[4];
    cout<< "SU2 TRANSPOSE IS USED"<<endl;
    for (int i = 0; i < 4; i++) { temp[i] = mat[i]; }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < i; j++)
        {
            mat[i*2+j] = temp[j*2+i];
            mat[j*2+i] = temp[i*2+j];
        }
    }
}

void SU2::conjugate()
{
    cout << "SU2 CONJUAGTE IS USED" << endl;
    for (int i = 0; i < 4; i++)
    {
        mat[i].conjugate();
    }
}

void SU2::copy(SU2 B)
{
    cout << "SU2 COPY IS USED" << endl;
    for (int i = 0; i < 4; i++)
    {
        mat[i] = B.mat[i];
    }
}

SU2 &SU2::operator+=(SU2 B)
{
    for (int i = 0; i < 4; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

SU2 &SU2::operator*=(SU2 B)
{
    /*
     * a*b = (a + bi)(c + id) = a*c + iad + ibc - bd;
     */
    complex temp[4];
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                temp[2*i+j] += mat[(2*i+k)]*B[(2*k+j)];
            }
//            temp[2*i+j] = mat[2*i]*B[j] + mat[2*i+1]*B[2+j];
        }
    }
//    temp[0] = mat[0]*B.mat[0] + mat[1]*B.mat[2];
//    temp[1] = mat[0]*B.mat[1] + mat[1]*B.mat[3];
//    temp[2] = mat[2]*B.mat[0] + mat[3]*B.mat[2];
//    temp[3] = mat[2]*B.mat[1] + mat[3]*B.mat[3];
    for (int i = 0; i < 4; i++) {
        mat[i] = temp[i];
    }
//    for (int i = 0; i < 2; i++) // ..this is faster?
//    {
//        for (int j = 0; j < 2; j++)
//        {
//            mat[(i*2+j)] = temp[(i*2+j)];
//        }
//    }
    return *this;
}

SU2 &SU2::operator*=(double b)
{
    for (int i = 0; i < 4; i++) {
        mat[i] *= b;
    }
    return *this;
}

SU2 &SU2::operator-=(SU2 B)
{
    for (int i = 0; i < 4; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

void SU2::zeros()
{
    for (int i = 0; i < 4; i++)
    {
        mat[i].setRe(0);
        mat[i].setIm(0);
    }
}

void SU2::identity()
{
    for (int i = 0; i < 4; i++)
    {
        mat[i].setIm(0);
    }
    mat[0].setRe(1);
    mat[1].setRe(0);
    mat[2].setRe(0);
    mat[3].setRe(1);
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
            cout << std::setw(15) << mat[2*i+j];
        }
        cout << endl;
    }
}
