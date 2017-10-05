#ifndef SU3_H
#define SU3_H

#include "complex.h"

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
    SU3 &operator/=(double a);
    void print();
    complex get(int i, int j) { return complex(mat[6*i + 2*j],mat[6*i + 2*j+1]); }

//    void copy(SU3 B);
    SU3 inv();
    void zeros();
    void identity();
    void transpose();
    SU3 conjugate();
    void copy(SU3 B);
    void setComplex(complex w, int i);

    // Complex number operations
    double norm(int i);
    double normSquared(int i);
    void c(int i);
    double re(int i) const { return mat[i]; } // But why constant?
    double im(int i) const { return mat[i+1]; }
    void setRe(int i, double re) { mat[i] = re; } // A bit redundant?
    void setIm(int i, double im) { mat[i+1] = im; }

//    complex(double real, double imag);
//    complex (const complex &b); // Copy constructor



//    // TEMP FOR PRINTING; MUST REMOVE TO STRIP DOWN LATER
//    friend std::ostream& operator<<(std::ostream& os, const complex& a); // Allows cout << myVector << endl;
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

#endif // SU3_H
