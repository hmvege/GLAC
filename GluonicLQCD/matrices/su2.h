#ifndef SU2_H
#define SU2_H

#include "complex.h"

class SU2
{
public:
    SU2();
    ~SU2();

//    complex *mat;
    complex mat[4];
    void print(); // TEMP, remove or comment out when program is done.
    complex get(int i, int j) { return mat[(2*i + j)]; }
    void copy(SU2 B);
    void zeros();
    void transpose();
    void conjugate();
//    SU2 inverse();
    complex &operator[](int i) { return mat[i]; }

    SU2 &operator+=(SU2 B);
    SU2 &operator-=(SU2 B);
    SU2 &operator*=(SU2 B);
    SU2 &operator*=(double B);
};

inline SU2 operator+(SU2 A, SU2 B)
{
    A += B;
    return A;
}

inline SU2 operator-(SU2 A, SU2 B)
{
    A -= B;
    return A;
}

inline SU2 operator*(SU2 A, SU2 B)
{
    A *= B;
    return A;
}

inline SU2 operator*(SU2 A, double b)
{
    A *= b;
    return A;
}


#endif // SU2_H
