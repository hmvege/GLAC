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

    double im;
    double re;
    complex &operator+=(complex b);
    complex &operator-=(complex b);
    complex &operator*=(complex b);
    complex &operator/=(complex b);
    complex &operator/=(double b);

    double norm();
    complex conjugate();

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
