#ifndef SU2_H
#define SU2_H

//#include <string>
#include <iostream>
#include <cmath>
#include "math/complex.h"

using std::cout;
using std::endl;

class SU2
{
public:
    SU2();
    ~SU2();

    double mat[8];

    void print();
    void zeros();
    void identity();
    void setComplex(complex w, int i);
    double normSquared(int i);

    SU2 transpose();
    SU2 conjugate();
    SU2 inv();

    // Getters, mostly using .mat[]
    double get(int i, int j, int k) { return mat[(4*i + 2*j) + k]; } // k is complex number of a position.
    double &operator[](int i) { return mat[i]; }

    // Operations
    SU2 &operator=(const SU2 &B);
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
