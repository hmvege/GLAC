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

//    // Copy constructor
//    SU3 (const SU3& other)
//    {
////        for (int i = 0; i < 18; i++) mat[i] = other.mat[i];
//        std::memcpy(&mat,&other.mat,144); // sizeof(double)*18 = 144
//    }

//    // Move constructor
//    SU3 (SU3 &&other) noexcept : mat{other.mat}
//    {
////        for (int i = 0; i < 18; i++) mat[i] = other.mat[i];
//    }

    // Copy assignement operator
    SU3 &operator =(const SU3& other)
    {
//        SU3 tmp(other);
//        *this = std::move(tmp);
//        return *this;

        for (int i = 0; i < 18; i++)
        {
            mat[i] = other.mat[i];
        }
        return *this;
    }

//    // Move assignement operator
//    SU3&operator= (SU3&& other) noexcept {
//        for (int i = 0; i < 18; i++)
//        {
//            mat[i] = other.mat[i];
//        }
//        return *this;
//    }

    // Element retrieval overloading
    double &operator[](int i) { return mat[i]; }

    // SU3 matrix operations overloading
    SU3 &operator+=(SU3 B);
    SU3 &operator-=(SU3 B);
    SU3 &operator*=(SU3 B);
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


#endif // SU3_H
