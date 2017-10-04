#ifndef _SU3_H
#define _SU3_H


class _su3
{
public:
    _su3();
    ~_su3();

    double mat[18];

    // Matrix specific functionss
    double &operator[](int i) { return mat[i]; }
    _su3 &operator =(const _su3 &B);
    _su3 &operator+=(_su3 B);
    _su3 &operator-=(_su3 B);
    _su3 &operator*=(_su3 B);
//    void print(); // TEMP, remove or comment out when program is done.
//    complex get(int i, int j) { return mat[(3*i + j)]; }

//    void copy(SU3 B);
    _su3 inv();
    void zeros();
    void identity();
    void transpose();
    void conjugate();
    void copy(_su3 B);

//    // Complex number operations
//    void c(int i); // Conjugate of complex number
    double re(int i);
    double im(int i);

//    complex(double real, double imag);

//    complex (const complex &b); // Copy constructor

//    double re(int i) const { return z[0]; } // But why constant?
//    double im(int i) const { return z[1]; }
//    void setRe(double re) { z[0] = re; }
//    void setIm(double im) { z[1] = im; }
////    void set(complex a);

//    complex &operator =(const complex& b);
//    complex &operator+=(complex b);
//    complex &operator-=(complex b);
//    complex &operator*=(complex b);
//    complex &operator*=(double b);
//    complex &operator/=(complex b);
//    complex &operator/=(double b);

//    double norm();
//    double normSquared();
//    complex conjugate();
//    complex c();

//    // TEMP FOR PRINTING; MUST REMOVE TO STRIP DOWN LATER
//    friend std::ostream& operator<<(std::ostream& os, const complex& a); // Allows cout << myVector << endl;
};

inline _su3 operator+(_su3 A, _su3 B)
{
    A += B;
    return A;
}

inline _su3 operator-(_su3 A, _su3 B)
{
    A -= B;
    return A;
}

inline _su3 operator*(_su3 A, _su3 B)
{
    A *= B;
    return A;
}

#endif // _SU3_H
