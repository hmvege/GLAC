#ifndef _SU2_H
#define _SU2_H

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>


class _SU2
{
public:
    _SU2();
    ~_SU2();

    double mat[8];

    void print();
    void copy(_SU2 B);
    void zeros();
    void identity();

    // Getters, mostly using .mat[]
    double get(int i, int j, int k) { return mat[(4*i + 2*j) + k]; } // k is complex number of a position.
    double &operator[](int i) { return mat[i]; }

    // Operations
    _SU2 &operator+=(_SU2 B);
    _SU2 &operator-=(_SU2 B);
    _SU2 &operator*=(_SU2 B);
    _SU2 &operator*=(double B);
};


inline _SU2 operator+(_SU2 A, _SU2 B)
{
A += B;
return A;
}

inline _SU2 operator-(_SU2 A, _SU2 B)
{
A -= B;
return A;
}

inline _SU2 operator*(_SU2 A, _SU2 B)
{
A *= B;
return A;
}

inline _SU2 operator*(_SU2 A, double b)
{
A *= b;
return A;
}

#endif // _SU2_H
