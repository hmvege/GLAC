#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include "math/complex.h"

template <class T>
class Lattice
{
private:
    std::vector<T> m_sites;
public:
    std::vector<int> m_dim; // Lattice dimensions
    int m_latticeSize;
    // Default contructor
    Lattice() {}
    Lattice(std::vector<int>latticeDimensions) {
        allocate(latticeDimensions);
    }

    // Destructor
    ~Lattice() {}

    // Copy constructor
    Lattice(const Lattice<T>& other) :
        m_dim(other.m_dim),
        m_latticeSize(other.m_latticeSize)
    {
        m_sites.resize(m_latticeSize);
        std::memcpy(&m_sites[0],&other.m_sites[0],sizeof(T)*m_latticeSize);
    }

    // Move constructor
    Lattice(Lattice<T> && other) noexcept :
        m_sites(other.m_sites),
        m_dim(other.m_dim),
        m_latticeSize(other.m_latticeSize)
    {
    }

    // Copy assignement operator
    Lattice &operator =(const Lattice& other) {
        Lattice tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    // Move assignement operator
    Lattice &operator= (Lattice<T>&& other) noexcept {
        m_latticeSize = other.m_latticeSize;
        m_dim  = other.m_dim;
        m_sites = other.m_sites;
        return *this;
    }

    void allocate(std::vector<int> dim);

    // Overloading lattice position getter
    T &operator[](int i) { return m_sites[i]; }
    // Overloading lattice operations
    Lattice<T> &operator+=(Lattice<T> B);
    Lattice<T> &operator-=(Lattice<T> B);
    Lattice<T> &operator*=(Lattice<T> B);
    // Operators using doubles, operations affect the whole of lattice
    Lattice<T> &operator*=(double B);
    Lattice<T> &operator/=(double B);
    // Operations involving complex operators
    Lattice<T> &operator+=(complex B);
    Lattice<T> &operator-=(complex B);
    Lattice<T> &operator*=(complex B);
    Lattice<T> &operator/=(complex B);

    // Lattice based operations SHOULD NOT BE CLASS MEMBERS, BUT HEADER FUNCTIONS!!!
    Lattice<double> realTrace();
    Lattice<double> imagTrace();
    Lattice<complex> trace();
    Lattice<T> sum();
    // Lattice value setters
    void identity();
    void zeros();
};

//////////////////////////////////////////
////// External operator overloading /////
//////////////////////////////////////////
// External lattice operator overloading
template <class T>
inline Lattice<T> operator+(Lattice<T> A, Lattice<T> B) {
    A += B;
    return A;
}

template <class T>
inline Lattice<T> operator-(Lattice<T> A, Lattice<T> B) {
    A -= B;
    return A;
}

template <class T>
inline Lattice<T> operator*(Lattice<T> A, Lattice<T> B) {
    A *= B;
    return A;
}

// External double operator overloading
template <class T>
inline Lattice<T> operator*(Lattice<T> A, double b) {
    A *= b;
    return A;
}

template <class T>
inline Lattice<T> operator/(Lattice<T> A, double b) {
    A /= b;
    return A;
}

//template <>
//inline Lattice<double> operator+(Lattice<double> A, b) {
//    A += b;
//}

// External complex operator overloading
template <class T>
inline Lattice<T> operator+(Lattice<T> A, complex b) {
    A += b;
    return A;
}

template <class T>
inline Lattice<T> operator-(Lattice<T> A, complex b) {
    A -= b;
    return A;
}

template <class T>
inline Lattice<T> operator*(Lattice<T> A, complex b) {
    A *= b;
    return A;
}

template <class T>
inline Lattice<T> operator/(Lattice<T> A, complex b) {
    A /= b;
    return A;
}

//////////////////////////////////////////
/////// Class operator overloading ///////
//////////////////////////////////////////
// Lattice operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator+=(Lattice<T> B) {
    printf("\nAdding!");
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] += B[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(Lattice<T> B) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] -= B[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(Lattice<T> B) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] *= B[iSite];
    }
    return *this;
}

// Double operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator*=(double b) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] *= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator/=(double b) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] /= b;
    }
    return *this;
}

// Complex operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator+=(complex b) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] += b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(complex b) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] -= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(complex b) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] *= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator/=(complex b) {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite] /= b;
    }
    return *this;
}

//////////////////////////////////////////
////////// Lattice operations ////////////
//////////////////////////////////////////
template <class T>
inline Lattice<double> Lattice<T>::realTrace()
{
    Lattice<double> tempTraceSum(m_dim);
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        tempTraceSum[iSite] = (m_sites[iSite][0] + m_sites[iSite][8] + m_sites[iSite][16]);
    }
    return tempTraceSum;
}

template <class T>
inline Lattice<double> Lattice<T>::imagTrace()
{
    Lattice<double> tempTraceSum(m_dim);
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        tempTraceSum[iSite] = (m_sites[iSite][1] + m_sites[iSite][9] + m_sites[iSite][17]);
    }
    return tempTraceSum;
}

template <class T>
inline Lattice<complex> Lattice<T>::trace()
{
    Lattice<complex> tempTraceSum(m_dim);
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        tempTraceSum[iSite] = complex(m_sites[iSite][0] + m_sites[iSite][8] + m_sites[iSite][16],
                                      m_sites[iSite][1] + m_sites[iSite][9] + m_sites[iSite][17]);
    }
    return tempTraceSum;
}

template <class T>
inline Lattice<T> Lattice<T>::sum()
{
    T latticeSum = 0;
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        latticeSum += m_sites[iSite];
    }
    return latticeSum;
}

//////////////////////////////////////////
///////// Lattice value setters //////////
//////////////////////////////////////////

template <class T>
inline void Lattice<T>::identity() {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite].identity();
    }
}

template <class T>
inline void Lattice<T>::zeros() {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite].zeros();
    }
}

// Allocates memory to the lattice. Has to be called every time(unless we are copying)
template <class T>
inline void Lattice<T>::allocate(std::vector<int> dim) {
    m_latticeSize = dim[0] * dim[1] * dim[2] * dim[3];
    m_dim = dim;
    m_sites.resize(m_latticeSize);
}

#endif // LATTICE_H
