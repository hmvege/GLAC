#ifndef SU3_H
#define SU3_H

#include "complex.h"

class SU3
{
public:
    SU3();
    ~SU3();
    SU3 &operator =(const SU3 &B);
    complex mat[9];
    void print(); // TEMP, remove or comment out when program is done.
    complex get(int i, int j) { return mat[(3*i + j)]; }

    void copy(SU3 B);
    void zeros();
    void transpose();
    void conjugate();
    SU3 inv();
    complex &operator[](int i) { return mat[i]; }

    SU3 &operator+=(SU3 B);
    SU3 &operator-=(SU3 B);
    SU3 &operator*=(SU3 B);
};

inline SU3 operator+(SU3 A, SU3 B)
{
    A += B;
    return A;
}

inline SU3 operator-(SU3 A, SU3 B)
{
    A -= B;
    return A;
}

inline SU3 operator*(SU3 A, SU3 B)
{
    A *= B;
    return A;
}

#endif // SU3_H
