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
    int m_N, m_Nx, m_Ny, m_Nz, m_Nt;
    // Default contructor
    Lattice() {}

    // Destructor
    ~Lattice() {}

    // Copy constructor
    Lattice(const Lattice<T>& other)
    {
//        printf("\nPRE: %d\n",m_N);
//        std::memcpy(&m_N,&other.m_N,sizeof(int));
        m_N = other.m_N;
        m_Nx = other.m_Nx;
        m_Ny = other.m_Ny;
        m_Nz = other.m_Nz;
        m_Nt = other.m_Nt;
//        allocate(m_N,m_Nx,m_Ny,m_Nz,m_Nt);
//        std::memcpy(&m_sites,&other.m_sites,sizeof(other.m_sites));
//        printf("\n%lu",sizeof(T));
        std::memcpy(&m_sites,&other.m_sites,sizeof(T)*m_N); // This should be 144*2, sizeof(T)*m_N
        printf("\nCopy constructor N = %d bytes = %lu\n",m_N,m_N*sizeof(T));
    }

    // Move constructor
    Lattice(Lattice<T> && other) noexcept :
        m_sites(other.m_sites),
        m_N(other.m_N),
        m_Nx(other.m_Nx),
        m_Ny(other.m_Ny),
        m_Nz(other.m_Nz),
        m_Nt(other.m_Nt)
    {
        printf("\nMove constructor\n");
    }

    // Copy assignement operator
    Lattice &operator =(const Lattice<T>& other) {
        printf("\nCopy assignement operator\n");
        Lattice tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    // Move assignement operator
    Lattice &operator= (Lattice<T>&& other) noexcept {
        m_N = other.m_N;
        m_Nx = other.m_Nx;
        m_Ny = other.m_Ny;
        m_Nz = other.m_Nz;
        m_Nt = other.m_Nt;
//        m_sites.resize(m_N);
        m_sites = other.m_sites;
//        other.m_sites.clear();
        printf("\nMove assignement operator\n");
        return *this;
    }

    void allocate(int N, int Nx, int Ny, int Nz, int Nt);

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
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] += B[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(Lattice<T> B) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] -= B[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(Lattice<T> B) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] *= B[iSite];
    }
    return *this;
}

// Double operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator*=(double b) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] *= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator/=(double b) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] /= b;
    }
    return *this;
}

// Complex operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator+=(complex b) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] += b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(complex b) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] -= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(complex b) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] *= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator/=(complex b) {
    for (int iSite = 0; iSite < m_N; iSite++) {
        m_sites[iSite] /= b;
    }
    return *this;
}


// Allocates memory to the lattice. Has to be called every time(unless we are copying)
template <class T>
inline void Lattice<T>::allocate(int N, int Nx, int Ny, int Nz, int Nt) {
    m_N = N;
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
    m_sites.resize(m_N);
}

#endif // LATTICE_H
