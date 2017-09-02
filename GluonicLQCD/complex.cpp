#include "complex.h"
#include <cmath>
using std::fabs;

/*
 * For storing complex numbers in the SU3 matrix class
 */

complex::complex() // Non-assigning constructor
{
//    re = 0;
//    im = 0;
    z = new double[2];
    z[0] = 0; // re
    z[1] = 0; // im
}

complex::~complex()
{
    delete [] z;
}
complex::complex(double real, double imag) // Value specific constructor

{
    z = new double[2];
    z[0] = real;
    z[1] = imag;
//    re = real;
//    im = imag;
}

complex::complex(const complex &b ) { // Copy constructor
    z = new double[2];
    z[0] = b.re();
    z[1] = b.im();
}

double complex::norm()
{
//    return re*re + im*im;
    return z[0]*z[0] + z[1]*z[1];
}

complex &complex::operator=(const complex &b)
{
//    if (this == &b)
//        return *this;

//    if (z) delete [] z; // Freeing up memory for new asignement
//    z = new double[2];

    z[0] = b.re();
    z[1] = b.im();
    return *this;
}


complex complex::conjugate()
{
//    im = -im;
    z[1] = -z[1];
    return *this;
}

complex complex::c()
{
    return complex(z[0],-z[1]);
//    return complex(re,-im);
}

complex &complex::operator+=(complex b)
{
//    re += b.re;
//    im += b.im;
    z[0] += b.z[0];
    z[1] += b.z[1];
    return *this;
}

complex &complex::operator*=(complex b)
{
    /*
     * a*b = (a + bi)(c + id) = a*c + iad + ibc - bd;
     */
//    double prev_re = re;
//    re = prev_re*b.re - im*b.im;
//    im = prev_re*b.im + im*b.re;
    double prev_re = z[0];
    z[0] = prev_re*b.re() - z[1]*b.im();
    z[1] = prev_re*b.im() + z[1]*b.re();
    return *this;
}

complex &complex::operator*=(double b)
{
    /*
     * a*b = (a + bi)(c + id) = a*c + iad + ibc - bd;
     */
//    double prev_re = re;
//    re = prev_re*b.re - im*b.im;
//    im = prev_re*b.im + im*b.re;
    z[0] *= b;
    z[1] *= b;
    return *this;
}

complex &complex::operator/=(complex b)
{
    /*
     * Dividing this/b
     */
//    double prev_re = re;
//    double divisor = b.re*b.re + b.im*b.im;
//    re = (re*b.re + im*b.im)/divisor;
//    im = (im*b.re - prev_re*b.im)/divisor;
    double prev_re = z[0];
    double divisor = b.re()*b.re() + b.im()*b.im();
    z[0] = (z[0]*b.re() + z[1]*b.im())/divisor;
    z[1] = (z[1]*b.re() - prev_re*b.im())/divisor;
    return *this;
}

complex &complex::operator/=(double b)
{
    /*
     * Dividing this/b
     */
    z[0] /= b;
    z[1] /= b;
    return *this;
}

complex &complex::operator-=(complex b)
{
    z[0] -= b.re();
    z[1] -= b.im();
    return *this;
}

// TEMP FOR PRINTING; MUST REMOVE TO STRIP DOWN LATER
std::ostream &operator<<(std::ostream &os, const complex &a)
{
    if (a.z[1] < 0) {
        os << a.z[0] << " - " << fabs(a.z[1]) << "i";
    }
    else {
        os << a.z[0] << " + " << a.z[1] << "i";
    }
    return os;
}
