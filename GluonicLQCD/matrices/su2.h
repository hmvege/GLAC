#ifndef SU2_H
#define SU2_H

#include "complex.h"

class SU2
{
public:
    SU2();
    ~SU2();

    complex mat[4];
    void print();
    void copy(SU2 B);
    void zeros();

    // Matrix operations
    void transpose();
    void conjugate();

    // Getters, mostly using .mat[]
    complex get(int i, int j) { return mat[(2*i + j)]; }
    complex &operator[](int i) { return mat[i]; }

    // Operations
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
