#ifndef SU3_H
#define SU3_H

#include "complex.h"

class SU3
{
//private:
//    static int matLength = 9;
public:
    SU3();
    ~SU3();
    complex mat[9];
    void print(); // TEMP
    complex get(int i, int j) { return mat[(3*i + j)]; }

    complex &operator[](int i) { return mat[i]; }

    SU3 &operator+=(SU3 B);
    SU3 &operator-=(SU3 B);
    SU3 &operator*=(SU3 B);
};

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

//inline int index(int i, int j)
//{
//    /*
//     * Function for contigious memory allocation.
//     */
//    return (3*i + j) ;
//}


#endif // SU3_H
