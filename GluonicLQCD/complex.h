#ifndef COMPLEX_H
#define COMPLEX_H

// TEMP FOR PRINTING; MUST REMOVE TO STRIP DOWN LATER
#include <string>
#include <iostream>

class complex
{
public:
    complex();
    ~complex();
    complex(double real, double imag);

    complex (const complex &b); // Copy constructor

    double *z;
//    double im;
//    double re;
    double re() const { return z[0]; } // But why constant?
    double im() const { return z[1]; }
    void setRe(double re) { z[0] = re; }
    void setIm(double im) { z[1] = im; }
//    void set(complex a);

    complex &operator =(const complex& b);
    complex &operator+=(complex b);
    complex &operator-=(complex b);
    complex &operator*=(complex b);
    complex &operator*=(double b);
    complex &operator/=(complex b);
    complex &operator/=(double b);

    double norm();
    complex conjugate();
    complex c();

    // TEMP FOR PRINTING; MUST REMOVE TO STRIP DOWN LATER
    friend std::ostream& operator<<(std::ostream& os, const complex& a); // Allows cout << myVector << endl;
};

inline complex operator+(complex a, complex b)
{
    a += b;
    return a;
}

inline complex operator-(complex a, complex b)
{
    a -= b;
    return a;
}

inline complex operator*(complex a, complex b)
{
    a *= b;
    return a;
}

inline complex operator*(complex a, double b)
{
    a *= b;
    return a;
}

inline complex operator/(complex a, complex b)
{
    a /= b;
    return a;
}

inline complex operator/(complex a, double b)
{
    a /= b;
    return a;
}


#endif // COMPLEX_H
