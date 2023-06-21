/*!
 * \class SU2
 *
 * \brief Class for holding \f$\mathrm{SU}(2)\f$ matrices.
 *
 * Operations have been overloaded for ease of use.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef SU2_H
#define SU2_H

#include "math/complex.h"

using std::cout;
using std::endl;

class SU2
{
public:
    SU2() {
        for (int i = 0; i < 8; i++)
        {
            mat[i] = 0;
        }
    }

    ~SU2() {}

    double mat[8];

    // Printers
    void print() const;
    friend std::ostream& operator<<(std::ostream& os, const SU2& mat);
    
    void zeros();
    void identity();
    void setComplex(complex w, int i);
    double normSquared(int i);

    SU2 transpose();
    SU2 conjugate();
    SU2 inv();

    // Getters, mostly using .mat[]
    inline double get(int i, int j, int k) { return mat[(4*i + 2*j) + k]; } // k is complex number of a position.
    inline double &operator[](int i) { return mat[i]; }

    // Comparison operators
    inline bool operator==(const SU2& other) const;
    inline bool operator!=(const SU2& other) const;

    // Operations
    SU2 &operator=(const SU2 &B);
    // Basic operations overloading with itself
    SU2 &operator+=(SU2 B);
    SU2 &operator-=(SU2 B);
    SU2 &operator*=(SU2 B);

    SU2 &operator*=(double B);
};

///////////////////////////////////////
//////// Operator overloading /////////
///////////////////////////////////////
inline SU2 operator+(SU2 A, SU2 B)
{
A += B;
return A;
}

inline SU2 operator-(SU2 A, SU2 B)
{
A -= B;
return A;
}

inline SU2 operator*(SU2 A, SU2 B)
{
A *= B;
return A;
}

inline SU2 operator*(SU2 A, double b)
{
A *= b;
return A;
}

///////////////////////////////
//// SU2 operator functions ///
///////////////////////////////

inline SU2 &SU2::operator=(const SU2 &B)
{
    /*
     * Copy assignement operator.
     */
    for (int i = 0; i < 8; i++) {
        mat[i] = B.mat[i];
    }
    return *this;
}

inline SU2 &SU2::operator+=(SU2 B)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] += B.mat[i];
    }
    return *this;
}

inline SU2 &SU2::operator-=(SU2 B)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] -= B.mat[i];
    }
    return *this;
}

inline SU2 &SU2::operator*=(SU2 B)
{
    double temp[8];

    temp[0] = mat[0]*B.mat[0] - mat[1]*B.mat[1] + mat[2]*B.mat[4] - mat[3]*B.mat[5];
    temp[1] = mat[0]*B.mat[1] + mat[1]*B.mat[0] + mat[2]*B.mat[5] + mat[3]*B.mat[4];
    temp[2] = mat[0]*B.mat[2] - mat[1]*B.mat[3] + mat[2]*B.mat[6] - mat[3]*B.mat[7];
    temp[3] = mat[0]*B.mat[3] + mat[1]*B.mat[2] + mat[2]*B.mat[7] + mat[3]*B.mat[6];
    temp[4] = mat[4]*B.mat[0] - mat[5]*B.mat[1] + mat[6]*B.mat[4] - mat[7]*B.mat[5];
    temp[5] = mat[4]*B.mat[1] + mat[5]*B.mat[0] + mat[6]*B.mat[5] + mat[7]*B.mat[4];
    temp[6] = mat[4]*B.mat[2] - mat[5]*B.mat[3] + mat[6]*B.mat[6] - mat[7]*B.mat[7];
    temp[7] = mat[4]*B.mat[3] + mat[5]*B.mat[2] + mat[6]*B.mat[7] + mat[7]*B.mat[6];

//    std::memcpy(mat, temp, 64);
    std::copy(temp, temp+8, mat);

    return *this;
}

inline SU2 &SU2::operator*=(double b)
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] *= b;
    }
    return *this;
}

inline bool SU2::operator==(const SU2 &other) const
{
    for (int i = 0; i < 8; i++)
    {
        if (mat[i] != other.mat[i])
            return false;
    }
    return true;
}

inline bool SU2::operator!=(const SU2 &other) const 
{
    return !(this->operator==(other));
}

///////////////////////////////
//// SU2 specific functions ///
///////////////////////////////

/*!
 * \brief SU2::inv
 * \return a copy of the inverse of itself.
 */
inline SU2 SU2::inv()
{
    SU2 R;
    R[0] = mat[0];
    R[1] = -mat[1];

    R[2] = mat[4];
    R[3] = -mat[5];

    R[4] = mat[2];
    R[5] = -mat[3];

    R[6] = mat[6];
    R[7] = -mat[7];

    return R;
}

///////////////////////////////
////// Matrix functions ///////
///////////////////////////////

inline void SU2::zeros()
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = 0;
    }
}

inline void SU2::identity()
{
    for (int i = 0; i < 8; i++)
    {
        mat[i] = 0;
    }
    mat[0] = 1;
    mat[6] = 1;
}

#endif // SU2_H
