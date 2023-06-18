/*!
 * \class complex
 *
 * \brief A complex number class, consisting of t
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>

class complex
{
public:
    complex();
    ~complex();
    complex(double real, double imag);

    complex (const complex &b); // Copy constructor

    double z[2];
    inline double re() const { return z[0]; }
    inline double im() const { return z[1]; }
    void setRe(double re) { z[0] = re; }
    void setIm(double im) { z[1] = im; }

    complex &operator =(const complex& b);

    // Equality operator
    inline bool operator==(const complex& a) const;
    inline bool operator!=(const complex& a) const;

    // Math operations
    complex &operator+=(complex b);
    complex &operator+=(double b);
    complex &operator-=(complex b);
    complex &operator-=(double b);
    complex &operator*=(complex b);
    complex &operator*=(double b);
    complex &operator/=(complex b);
    complex &operator/=(double b);
    complex operator-() const;

    double norm();
    double normSquared();
    complex conjugate();
    complex c();
    complex zeros();

    friend std::ostream& operator<<(std::ostream& os, const complex& a);
};

inline bool complex::operator==(const complex &a) const
{
    if (re() != a.re() || im() != a.im())
        return false;
    
    return true;
}

inline bool complex::operator!=(const complex &other) const 
{
    return !(this->operator==(other));
}

inline complex operator+(complex a, complex b)
{
    a += b;
    return a;
}

inline complex operator+(complex a, double b)
{
    a += b;
    return a;
}

inline complex operator-(complex a, complex b)
{
    a -= b;
    return a;
}

inline complex operator-(complex a, double b)
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
