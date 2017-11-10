#ifndef SU3_H
#define SU3_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include "math/complex.h"

using std::cout;
using std::endl;

class SU3
{
public:
    SU3();
    ~SU3();

    double mat[18];

    // Matrix specific functionss
    double &operator[](int i) { return mat[i]; }
    SU3 &operator =(const SU3 &B);
    SU3 &operator+=(SU3 B);
    SU3 &operator-=(SU3 B);
    SU3 &operator*=(SU3 B);
    SU3 &operator*=(complex z);
    SU3 &operator/=(double a);
    SU3 &operator*=(double a);
    SU3 &operator-=(double a);
    void print();
    void printMachine();
    complex get(int i, int j) { return complex(mat[6*i + 2*j],mat[6*i + 2*j+1]); }

    complex trace();
    SU3 makeHermitian();
    SU3 makeAntiHermitian();

    SU3 getIm();
    SU3 getRe();

    SU3 inv();
    void zeros();
    void identity();
    SU3 transpose();
    SU3 conjugate();
    void copy(SU3 B);
    void setComplex(complex w, int i);

    // Complex number operations
    double norm(int i);
    double normSquared(int i);
    complex c(int i);
    double re(int i) const { return mat[i]; } // But why constant?
    double im(int i) const { return mat[i+1]; }
    void setRe(int i, double re) { mat[i] = re; } // A bit redundant?
    void setIm(int i, double im) { mat[i+1] = im; }

    friend std::ostream& operator<<(std::ostream& os, const SU3& A); // Allows cout << myVector << endl;
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

inline SU3 operator/(SU3 A, double a)
{
    A /= a;
    return A;
}

inline SU3 operator*(SU3 A, double a)
{
    A *= a;
    return A;
}

inline SU3 operator-(SU3 A, double a)
{
    A -= a;
    return A;
}

inline SU3 operator*(SU3 A, complex z)
{
    A *= z;
    return A;
}

#endif // SU3_H
