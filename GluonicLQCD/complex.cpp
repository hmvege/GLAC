#include "complex.h"

/*
 * For storing complex numbers in the SU3 matrix class
 */

complex::complex()
{

}

complex::~complex()
{

}

complex::complex(double real, double imag)
{
    re = real;
    im = imag;
}

complex complex::conjugate()
{
    im = -im;
    return *this;
}

complex &complex::operator+=(complex b)
{
    re = re + b.re;
    im = im + b.im;
    return *this;
}

complex &complex::operator*=(complex b)
{
    /*
     * a*b = (a + bi)(c + id) = a*c + iad + ibc - bd;
     */
    double prev_re = re;
    re = prev_re*b.re - im*b.im;
    im = prev_re*b.im + im*b.re;
    return *this;
}

complex &complex::operator/=(complex b)
{
    /*
     * Dividing this/b
     */
    double prev_re = re;
    double divisor = b.re*b.re + b.im*b.im;
    re = (re*b.re + im*b.im)/divisor;
    im = (im*b.re - prev_re*b.im)/divisor;
    return *this;
}

complex &complex::operator-=(complex b)
{
    re = re - b.re;
    im = im - b.im;
    return *this;
}

// TEMP FOR PRINTING; MUST REMOVE TO STRIP DOWN LATER
std::ostream &operator<<(std::ostream &os, const complex &a)
{
    os << a.re << " + i" << a.im;
    return os;
}
