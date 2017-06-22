#include "complex.h"
#include <cmath>
using std::fabs;

/*
 * For storing complex numbers in the SU3 matrix class
 */

complex::complex()
{
    re = 0;
    im = 0;
}

complex::~complex()
{

}

double complex::norm()
{
    return re*re + im*im;
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

complex complex::c()
{
    return complex(re,-im);
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

complex &complex::operator/=(double b)
{
    /*
     * Dividing this/b
     */
    re /= b;
    im /= b;
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
    if (a.im < 0) {
        os << a.re << " - " << fabs(a.im) << "i";
    }
    else {
        os << a.re << " + " << a.im << "i";
    }
    return os;
}
