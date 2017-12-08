#ifndef SU3_H
#define SU3_H

#include "math/complex.h"

class SU3
{
public:
    // Static matrix allocation
    double mat[18];

    // Default contructor
    SU3() {}

    // Destructor
    ~SU3() {}

    // Copy assignement operator
    SU3 &operator =(const SU3& other)
    {
        for (int i = 0; i < 18; i++)
        {
            mat[i] = other.mat[i];
        }
        return *this;
    }

    // Element retrieval overloading
    double &operator[](int i) { return mat[i]; }

    // Overloading the setting operator
    SU3 &operator =(const double& other);

    // SU3 matrix operations overloading
    SU3 &operator+=(const SU3 &B);
    SU3 &operator+=(SU3 &&B);
    SU3 &operator-=(const SU3&B);
    SU3 &operator-=(SU3 && B);
    SU3 &operator*=(const SU3 &B);
    SU3 &operator*=(SU3 && B);
    // Complex matrix operations overloading
    SU3 &operator+=(complex z);
    SU3 &operator-=(complex z);
    SU3 &operator*=(complex z);
    // Double matrix operations overloading
    SU3 &operator/=(double a);
    SU3 &operator*=(double a);

    // Printing functions
    void print();

    // Complex element getter
    complex get(int i, int j) { return complex(mat[6*i + 2*j],mat[6*i + 2*j+1]); }

    // Complex trace, implement a real trace as well?
    complex trace();

    // Functions for making matrix hermitian/anti-hermitian
    SU3 makeHermitian();
    SU3 makeAntiHermitian();

    // Returns a matrix of only real/imaginary elements
    SU3 getIm();
    SU3 getRe();

    // Returns the inverse of the matrix(the conjugate transpose)
    SU3 inv();

    // Sets the matrix to zero
    void zeros();

    // Sets the matrix to identity
    void identity();

    // Transposes the matrix
    SU3 transpose();

    // Complex conjugates the matrix
    SU3 conjugate();

    // Sets position i (contigious allocation) as a complex number.
    void setComplex(complex w, int i);

    // Complex number operations
    double norm(int i);
    double normSquared(int i);
    complex c(int i);
    double re(int i) const { return mat[i]; } // But why constant?
    double im(int i) const { return mat[i+1]; }
    void setRe(int i, double re) { mat[i] = re; } // A bit redundant?
    void setIm(int i, double im) { mat[i+1] = im; }
};

///////////////////////////////////////
//////// Operator overloading /////////
///////////////////////////////////////
// SU3 operator overloading
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

// Double operator overloading
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

// Complex operator overloading
inline SU3 operator+(SU3 A, complex z)
{
    A += z;
    return A;
}

inline SU3 operator-(SU3 A, complex z)
{
    A -= z;
    return A;
}

inline SU3 operator*(SU3 A, complex z)
{
    A *= z;
    return A;
}

///////////////////////////////
//// Overloading operators ////
///////////////////////////////

inline SU3 &SU3::operator+=(const SU3& B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

inline SU3 &SU3::operator+=(SU3&& B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

inline SU3 &SU3::operator-=(const SU3 &B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

inline SU3 &SU3::operator-=(SU3 &&B)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

inline SU3 &SU3::operator*=(const SU3 &B)
{
    /*
     * ab = (a + bi)(c + id) = ac + iad + ibc - bd = ac - bd + i(ad + bc);
     * 0 1 2   11 12 13   00 01 02
     * 3 4 5 = 21 22 23 = 10 11 12
     * 6 7 8   31 32 33   20 21 22
     */
    double temp[18];

    temp[0] = mat[0]*B.mat[0] - mat[1]*B.mat[1] + mat[2]*B.mat[6] - mat[3]*B.mat[7] + mat[4]*B.mat[12] - mat[5]*B.mat[13];
    temp[1] = mat[0]*B.mat[1] + mat[1]*B.mat[0] + mat[2]*B.mat[7] + mat[3]*B.mat[6] + mat[4]*B.mat[13] + mat[5]*B.mat[12];
    temp[2] = mat[0]*B.mat[2] - mat[1]*B.mat[3] + mat[2]*B.mat[8] - mat[3]*B.mat[9] + mat[4]*B.mat[14] - mat[5]*B.mat[15];
    temp[3] = mat[0]*B.mat[3] + mat[1]*B.mat[2] + mat[2]*B.mat[9] + mat[3]*B.mat[8] + mat[4]*B.mat[15] + mat[5]*B.mat[14];
    temp[4] = mat[0]*B.mat[4] - mat[1]*B.mat[5] + mat[2]*B.mat[10] - mat[3]*B.mat[11] + mat[4]*B.mat[16] - mat[5]*B.mat[17];
    temp[5] = mat[0]*B.mat[5] + mat[1]*B.mat[4] + mat[2]*B.mat[11] + mat[3]*B.mat[10] + mat[4]*B.mat[17] + mat[5]*B.mat[16];
    temp[6] = mat[6]*B.mat[0] - mat[7]*B.mat[1] + mat[8]*B.mat[6] - mat[9]*B.mat[7] + mat[10]*B.mat[12] - mat[11]*B.mat[13];
    temp[7] = mat[6]*B.mat[1] + mat[7]*B.mat[0] + mat[8]*B.mat[7] + mat[9]*B.mat[6] + mat[10]*B.mat[13] + mat[11]*B.mat[12];
    temp[8] = mat[6]*B.mat[2] - mat[7]*B.mat[3] + mat[8]*B.mat[8] - mat[9]*B.mat[9] + mat[10]*B.mat[14] - mat[11]*B.mat[15];
    temp[9] = mat[6]*B.mat[3] + mat[7]*B.mat[2] + mat[8]*B.mat[9] + mat[9]*B.mat[8] + mat[10]*B.mat[15] + mat[11]*B.mat[14];
    temp[10] = mat[6]*B.mat[4] - mat[7]*B.mat[5] + mat[8]*B.mat[10] - mat[9]*B.mat[11] + mat[10]*B.mat[16] - mat[11]*B.mat[17];
    temp[11] = mat[6]*B.mat[5] + mat[7]*B.mat[4] + mat[8]*B.mat[11] + mat[9]*B.mat[10] + mat[10]*B.mat[17] + mat[11]*B.mat[16];
    temp[12] = mat[12]*B.mat[0] - mat[13]*B.mat[1] + mat[14]*B.mat[6] - mat[15]*B.mat[7] + mat[16]*B.mat[12] - mat[17]*B.mat[13];
    temp[13] = mat[12]*B.mat[1] + mat[13]*B.mat[0] + mat[14]*B.mat[7] + mat[15]*B.mat[6] + mat[16]*B.mat[13] + mat[17]*B.mat[12];
    temp[14] = mat[12]*B.mat[2] - mat[13]*B.mat[3] + mat[14]*B.mat[8] - mat[15]*B.mat[9] + mat[16]*B.mat[14] - mat[17]*B.mat[15];
    temp[15] = mat[12]*B.mat[3] + mat[13]*B.mat[2] + mat[14]*B.mat[9] + mat[15]*B.mat[8] + mat[16]*B.mat[15] + mat[17]*B.mat[14];
    temp[16] = mat[12]*B.mat[4] - mat[13]*B.mat[5] + mat[14]*B.mat[10] - mat[15]*B.mat[11] + mat[16]*B.mat[16] - mat[17]*B.mat[17];
    temp[17] = mat[12]*B.mat[5] + mat[13]*B.mat[4] + mat[14]*B.mat[11] + mat[15]*B.mat[10] + mat[16]*B.mat[17] + mat[17]*B.mat[16];

    for (int i = 0; i < 18; i++)
    {
        mat[i] = temp[i];
    }
    return *this;
}

inline SU3 &SU3::operator*=(SU3&& B)
{
    /*
     * ab = (a + bi)(c + id) = ac + iad + ibc - bd = ac - bd + i(ad + bc);
     * 0 1 2   11 12 13   00 01 02
     * 3 4 5 = 21 22 23 = 10 11 12
     * 6 7 8   31 32 33   20 21 22
     */
    double temp[18];

    temp[0] = mat[0]*B.mat[0] - mat[1]*B.mat[1] + mat[2]*B.mat[6] - mat[3]*B.mat[7] + mat[4]*B.mat[12] - mat[5]*B.mat[13];
    temp[1] = mat[0]*B.mat[1] + mat[1]*B.mat[0] + mat[2]*B.mat[7] + mat[3]*B.mat[6] + mat[4]*B.mat[13] + mat[5]*B.mat[12];
    temp[2] = mat[0]*B.mat[2] - mat[1]*B.mat[3] + mat[2]*B.mat[8] - mat[3]*B.mat[9] + mat[4]*B.mat[14] - mat[5]*B.mat[15];
    temp[3] = mat[0]*B.mat[3] + mat[1]*B.mat[2] + mat[2]*B.mat[9] + mat[3]*B.mat[8] + mat[4]*B.mat[15] + mat[5]*B.mat[14];
    temp[4] = mat[0]*B.mat[4] - mat[1]*B.mat[5] + mat[2]*B.mat[10] - mat[3]*B.mat[11] + mat[4]*B.mat[16] - mat[5]*B.mat[17];
    temp[5] = mat[0]*B.mat[5] + mat[1]*B.mat[4] + mat[2]*B.mat[11] + mat[3]*B.mat[10] + mat[4]*B.mat[17] + mat[5]*B.mat[16];
    temp[6] = mat[6]*B.mat[0] - mat[7]*B.mat[1] + mat[8]*B.mat[6] - mat[9]*B.mat[7] + mat[10]*B.mat[12] - mat[11]*B.mat[13];
    temp[7] = mat[6]*B.mat[1] + mat[7]*B.mat[0] + mat[8]*B.mat[7] + mat[9]*B.mat[6] + mat[10]*B.mat[13] + mat[11]*B.mat[12];
    temp[8] = mat[6]*B.mat[2] - mat[7]*B.mat[3] + mat[8]*B.mat[8] - mat[9]*B.mat[9] + mat[10]*B.mat[14] - mat[11]*B.mat[15];
    temp[9] = mat[6]*B.mat[3] + mat[7]*B.mat[2] + mat[8]*B.mat[9] + mat[9]*B.mat[8] + mat[10]*B.mat[15] + mat[11]*B.mat[14];
    temp[10] = mat[6]*B.mat[4] - mat[7]*B.mat[5] + mat[8]*B.mat[10] - mat[9]*B.mat[11] + mat[10]*B.mat[16] - mat[11]*B.mat[17];
    temp[11] = mat[6]*B.mat[5] + mat[7]*B.mat[4] + mat[8]*B.mat[11] + mat[9]*B.mat[10] + mat[10]*B.mat[17] + mat[11]*B.mat[16];
    temp[12] = mat[12]*B.mat[0] - mat[13]*B.mat[1] + mat[14]*B.mat[6] - mat[15]*B.mat[7] + mat[16]*B.mat[12] - mat[17]*B.mat[13];
    temp[13] = mat[12]*B.mat[1] + mat[13]*B.mat[0] + mat[14]*B.mat[7] + mat[15]*B.mat[6] + mat[16]*B.mat[13] + mat[17]*B.mat[12];
    temp[14] = mat[12]*B.mat[2] - mat[13]*B.mat[3] + mat[14]*B.mat[8] - mat[15]*B.mat[9] + mat[16]*B.mat[14] - mat[17]*B.mat[15];
    temp[15] = mat[12]*B.mat[3] + mat[13]*B.mat[2] + mat[14]*B.mat[9] + mat[15]*B.mat[8] + mat[16]*B.mat[15] + mat[17]*B.mat[14];
    temp[16] = mat[12]*B.mat[4] - mat[13]*B.mat[5] + mat[14]*B.mat[10] - mat[15]*B.mat[11] + mat[16]*B.mat[16] - mat[17]*B.mat[17];
    temp[17] = mat[12]*B.mat[5] + mat[13]*B.mat[4] + mat[14]*B.mat[11] + mat[15]*B.mat[10] + mat[16]*B.mat[17] + mat[17]*B.mat[16];

    for (int i = 0; i < 18; i++)
    {
        mat[i] = temp[i];
    }
    return *this;
}

inline SU3 &SU3::operator/=(double a)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] /= a;
    }
    return *this;
}

inline SU3 &SU3::operator*=(double a)
{
    for (int i = 0; i < 18; i++)
    {
        mat[i] *= a;
    }
    return *this;
}

inline SU3 &SU3::operator+=(complex z)
{
    for (int i = 0; i < 18; i+=2)
    {
        mat[i] += z.z[0];
        mat[i+1] += z.z[1];
    }
    return *this;
}

inline SU3 &SU3::operator-=(complex z)
{
    for (int i = 0; i < 18; i+=2)
    {
        mat[i] -= z.z[0];
        mat[i+1] -= z.z[1];
    }
    return *this;
}

inline SU3 &SU3::operator*=(complex z)
{
    double temp;
    for (int i = 0; i < 18; i+=2)
    {
        temp = mat[i];
        mat[i] = mat[i]*z.z[0] - mat[i+1]*z.z[1];
        mat[i+1] = temp*z.z[1] + mat[i+1]*z.z[0];
    }
    return *this;
}

inline SU3 &SU3::operator =(const double& other)
{
    for (int i = 0; i < 18; i++) {
        mat[i] = other;
    }
    return *this;
}

#endif // SU3_H
