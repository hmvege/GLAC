/*!
 * \class Lattice
 *
 * \brief A method for holding a lattice of type T.
 *
 * The lattice class is implemented such that is uses move semantics, i.e. rule of 5.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <mpi.h>
#include "parallelization/index.h"
#include "parallelization/neighbours.h"
#include "parallelization/parallelparameters.h"
#include "functions.h"

template <class T>
class Lattice
{
public:
    std::vector<T> m_sites;
    std::vector<unsigned int> m_dim; // Lattice dimensions
    unsigned long m_latticeSize;

    // Default contructors
    Lattice() {}
    Lattice(std::vector<unsigned int> latticeDimensions)
    {
        allocate(latticeDimensions);
    }

    // Destructor
    ~Lattice() {}

    // Copy constructor
    /*!
     * \brief Lattice copy constructor
     * \param other
     */
    Lattice(const Lattice<T> &other) :
        m_dim(other.m_dim),
        m_latticeSize(other.m_latticeSize)
    {
        m_sites = other.m_sites;
    }

    // Move constructor
    /*!
     * \brief Lattice move constructor
     * \param other
     */
    Lattice(Lattice<T> &&other) noexcept :
        m_sites(std::move(other.m_sites)),
        m_dim(other.m_dim),
        m_latticeSize(other.m_latticeSize)
    {
    }

    // Copy assignement operator
    /*!
     * \brief operator = Copy assignement operator
     * \param other
     * \return
     */
    Lattice &operator =(const Lattice& other)
    {
        Lattice tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    // Move assignement operator
    /*!
     * \brief operator = Move assignement operator
     * \param other
     * \return
     */
    Lattice &operator= (Lattice<T> &&other) noexcept
    {
        m_latticeSize = other.m_latticeSize;
        m_dim  = other.m_dim;
        m_sites = std::move(other.m_sites);
        return *this;
    }

    /*!
     * \brief allocate method for allocating memory.
     * \param dim
     */
    void allocate(const std::vector<unsigned int> &dim);

    // Overloading lattice position getter
    /*!
     * \brief operator [] access operator that has no boundary checking.
     * \param i
     * \return
     */
    T &operator[](const unsigned long i)
    {
        return m_sites[i];
    }

    /*!
     * \brief operator () access operator that has boundary checking.
     * \param i
     * \return
     */
    T &operator()(const unsigned long i)
    {
        if (i >= m_latticeSize)
        {
            printf("OUT OF BUNDS WITH INDEX %lu of %lu", i, m_latticeSize);
            exit(1);
        }
        return m_sites[i];
    }

    // Overloading lattice operations
    Lattice<T> &operator+=(const Lattice<T>& B);
    Lattice<T> &operator+=(Lattice<T>&& B);
    Lattice<T> &operator-=(const Lattice<T>& B);
    Lattice<T> &operator-=(Lattice<T>&& B);
    Lattice<T> &operator*=(const Lattice<T>& B);
    Lattice<T> &operator*=(Lattice<T>&& B);
    // Operators using doubles, operations affect the whole of lattice
    Lattice<T> &operator*=(const double B);
    Lattice<T> &operator/=(const double B);
    // Operations involving complex operators
    Lattice<T> &operator+=(const complex B);
    Lattice<T> &operator-=(const complex B);
    Lattice<T> &operator*=(const complex B);
    Lattice<T> &operator/=(const complex B);

    /*!
     * \brief identity method for setting the SU3 matrices of the lattice to identity.
     */
    void identity();

    /*!
     * \brief zeros sets all of the SU3 elements in the lattice to zero.
     */
    void zeros();

    /*!
     * \brief makeHermitian converts the SU3 matrices of the lattice to hermitian matrices.
     */
    void makeHermitian();

    /*!
     * \brief makeHermitian converts the SU3 matrices of the lattice to anti-hermitian matrices.
     */
    void makeAntiHermitian();

    /*!
     * \brief copy method for copying content of another lattice B onto itself.
     * \param B
     * \return itself as copied from B.
     */
    Lattice<T> copy(const Lattice<T> &B);

    // Make private?
//    static std::vector<SU3> sendCubeX;
//    static std::vector<SU3> recvCubeX;
//    static std::vector<SU3> sendCubeY;
//    static std::vector<SU3> recvCubeY;
//    static std::vector<SU3> sendCubeZ;
//    static std::vector<SU3> recvCubeZ;
//    static std::vector<SU3> sendCubeT;
//    static std::vector<SU3> recvCubeT;

//    void initializeShiftVariables() {
//        sendCubeX.resize();
//    }
};

//void initializseLatticeSharing(std::vector<unsigned int> dim) {
//    sendCubeX.resize(dim[1]*dim[2]*dim[3]);
//    recvCubeX.resize(dim[1]*dim[2]*dim[3]);
//    sendCubeY.resize(dim[0]*dim[2]*dim[3]);
//    recvCubeY.resize(dim[0]*dim[2]*dim[3]);
//    sendCubeZ.resize(dim[0]*dim[1]*dim[3]);
//    recvCubeZ.resize(dim[0]*dim[1]*dim[3]);
//    sendCubeT.resize(dim[0]*dim[1]*dim[2]);
//    recvCubeT.resize(dim[0]*dim[1]*dim[2]);
//}

//std::vector<SU3> Lattice::sendCubeX;
//std::vector<SU3> Lattice<SU3>::recvCubeX;
//std::vector<SU3> sendCubeY;
//std::vector<SU3> recvCubeY;
//std::vector<SU3> sendCubeZ;
//std::vector<SU3> recvCubeZ;
//std::vector<SU3> sendCubeT;
//std::vector<SU3> recvCubeT;

//////////////////////////////////////////
////// External operator overloading /////
//////////////////////////////////////////
// External lattice operator overloading
template <class T>
inline Lattice<T> operator+(Lattice<T> A, Lattice<T>& B)
{
    A += B;
    return A;
}

template <class T>
inline Lattice<T> operator+(Lattice<T> A, Lattice<T>&& B)
{
    A += B;
    return A;
}

template <class T>
inline Lattice<T> operator-(Lattice<T> A, Lattice<T>& B)
{
    A -= B;
    return A;
}

template <class T>
inline Lattice<T> operator-(Lattice<T> A, Lattice<T>&& B)
{
    A -= B;
    return A;
}

template <class T>
inline Lattice<T> operator*(Lattice<T> A, Lattice<T>& B)
{
    A *= B;
    return A;
}

template <class T>
inline Lattice<T> operator*(Lattice<T> A, Lattice<T>&& B)
{
    A *= B;
    return A;
}

// External double operator overloading
template <class T>
inline Lattice<T> operator*(Lattice<T> A, double b)
{
    A *= b;
    return A;
}

template <class T>
inline Lattice<T> operator/(Lattice<T> A, double b)
{
    A /= b;
    return A;
}

// External complex operator overloading
template <class T>
inline Lattice<T> operator+(Lattice<T> A, complex b)
{
    A += b;
    return A;
}

template <class T>
inline Lattice<T> operator-(Lattice<T> A, complex b)
{
    A -= b;
    return A;
}

template <class T>
inline Lattice<T> operator*(Lattice<T> A, complex b)
{
    A *= b;
    return A;
}

template <class T>
inline Lattice<T> operator/(Lattice<T> A, complex b)
{
    A /= b;
    return A;
}

//////////////////////////////////////////
/////// Class operator overloading ///////
//////////////////////////////////////////
// Lattice operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator+=(const Lattice<T>& B)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] += B.m_sites[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator+=(Lattice<T>&& B)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] += B.m_sites[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(const Lattice<T>& B)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] -= B.m_sites[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(Lattice<T>&& B)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] -= B.m_sites[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(const Lattice<T>& B)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] *= B.m_sites[iSite];
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(Lattice<T>&& B)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] *= B.m_sites[iSite];
    }
    return *this;
}

// Double operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator*=(const double b)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] *= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator/=(const double b)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] /= b;
    }
    return *this;
}

// Complex operator overloading
template <class T>
inline Lattice<T> &Lattice<T>::operator+=(const complex b)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] += b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator-=(const complex b)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] -= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator*=(const complex b)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] *= b;
    }
    return *this;
}

template <class T>
inline Lattice<T> &Lattice<T>::operator/=(const complex b)
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] /= b;
    }
    return *this;
}

//////////////////////////////////////////
///////// Lattice value setters //////////
//////////////////////////////////////////
// Value setters of the lattice
template <class T>
inline void Lattice<T>::identity()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite].identity();
    }
}

template <>
inline void Lattice<complex>::identity()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite].z[0] = 1.0;
        m_sites[iSite].z[1] = 0.0;
    }
}

template <>
inline void Lattice<double>::identity()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] = 1.0;
    }
}

template <class T>
inline void Lattice<T>::zeros()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite].zeros();
    }
}

template <>
inline void Lattice<complex>::zeros()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite].z[0] = 0;
        m_sites[iSite].z[1] = 0;
    }
}

template <>
inline void Lattice<double>::zeros()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] = 0;
    }
}

template <>
inline void Lattice<SU3>::makeHermitian()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite].makeHermitian();
    }
}

template <>
inline void Lattice<SU3>::makeAntiHermitian()
{
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite].makeAntiHermitian();
    }
}

// Allocates memory to the lattice. Has to be called every time(unless we are copying)
template <class T>
inline void Lattice<T>::allocate(const std::vector<unsigned int> &dim)
{
    m_latticeSize = dim[0] * dim[1] * dim[2] * dim[3];
    m_dim = dim;
    m_sites.resize(m_latticeSize);
}

template <class T> // TEMP TEST!!
inline Lattice<T> Lattice<T>::copy(const Lattice<T> &B)
{
    m_dim = B.m_dim;
    m_latticeSize = B.m_latticeSize;
    m_sites.resize(m_latticeSize);
    for (unsigned long int iSite = 0; iSite < m_latticeSize; iSite++)
    {
        m_sites[iSite] = B.m_sites[iSite];
    }
    return *this;
}

//////////////////////////////////////////
////////// Lattice functions /////////////
//////////////////////////////////////////
template <class SU3>
inline Lattice<double> realTrace(const Lattice<SU3> &L)
{
    Lattice<double> tempTraceSum(L.m_dim);
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        tempTraceSum[iSite] = (L.m_sites[iSite].mat[0] + L.m_sites[iSite].mat[8] + L.m_sites[iSite].mat[16]);
    }
    return tempTraceSum;
}

template <class SU3>
inline Lattice<double> imagTrace(const Lattice<SU3> &L)
{
    Lattice<double> tempTraceSum(L.m_dim);
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        tempTraceSum[iSite] = (L.m_sites[iSite].mat[1] + L.m_sites[iSite].mat[9] + L.m_sites[iSite].mat[17]);
    }
    return tempTraceSum;
}

template <class SU3>
inline Lattice<complex> trace(const Lattice<SU3> &L)
{
    Lattice<complex> tempTraceSum(L.m_dim);
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        tempTraceSum[iSite] = complex(L.m_sites[iSite].mat[0] + L.m_sites[iSite].mat[8] + L.m_sites[iSite].mat[16],
                                      L.m_sites[iSite].mat[1] + L.m_sites[iSite].mat[9] + L.m_sites[iSite].mat[17]);
    }
    return tempTraceSum;
}

template <class SU3>
inline Lattice<SU3> subtractImag(Lattice <SU3> &L, const Lattice <double> &other)
{
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        L[iSite][1] -= other.m_sites[iSite];
        L[iSite][9] -= other.m_sites[iSite];
        L[iSite][17] -= other.m_sites[iSite];
    }
    return std::move(L);
}

template <class SU3>
inline Lattice<SU3> subtractImag(Lattice <SU3> &&L, Lattice <double> &&other)
{
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        L[iSite][1] -= other[iSite];
        L[iSite][9] -= other[iSite];
        L[iSite][17] -= other[iSite];
    }
    return std::move(L);
}

template <class SU3>
inline Lattice<SU3> subtractReal(Lattice <SU3> &L, const Lattice <double> &other)
{
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        L[iSite][0] -= other.m_sites[iSite];
        L[iSite][8] -= other.m_sites[iSite];
        L[iSite][16] -= other.m_sites[iSite];
    }
    return std::move(L);
}

template <class SU3>
inline Lattice<SU3> subtractReal(Lattice <SU3> &&L, Lattice <double> &&other)
{
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        L[iSite][0] -= other[iSite];
        L[iSite][8] -= other[iSite];
        L[iSite][16] -= other[iSite];
    }
    return std::move(L);
}

template <class T>
inline T sum(const Lattice<T> &L)
{
    T latticeSum;
    latticeSum = 0.0;
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        latticeSum += L.m_sites[iSite];
    }
    return latticeSum;
}

//template <>
inline std::vector<double> sumSpatial(const Lattice<double> &L)
{
    /*
     * For summing the topological charge xyz into the time axis.
     */

    // Creates empty vector for time axis points
    std::vector<double> latticeSpatialSum(L.m_dim[3], 0);

    // Sums the xyz directions into the time axis
    for (unsigned int it = 0; it < L.m_dim[3]; it++)
    {
        for (unsigned int iy = 0; iy < L.m_dim[1]; iy++)
        {
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++)
            {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++)
                {
                    latticeSpatialSum[it] += L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                }
            }
        }
    }

    return latticeSpatialSum;
}

inline Lattice<double> realTraceMultiplication(const Lattice<SU3> &L1, const Lattice<SU3> &L2)
{
    /*
     * Function for multiplying lattice and taking the real trace without summing the matrix
     */
    unsigned long site = 0;
    Lattice<double> latticeSum;
    latticeSum.allocate(L1.m_dim);
    latticeSum.zeros();

    for (unsigned int ix = 0; ix < L1.m_dim[0]; ix++)
    {
        for (unsigned int iy = 0; iy < L1.m_dim[1]; iy++)
        {
            for (unsigned int iz = 0; iz < L1.m_dim[2]; iz++)
            {
                for (unsigned int it = 0; it < L1.m_dim[3]; it++)
                {
                    site = Parallel::Index::getIndex(ix,iy,iz,it);
                    latticeSum[site] += traceRealMultiplication(L1.m_sites[site], L2.m_sites[site]);
                }
            }
        }
    }

    return latticeSum;
}


inline Lattice<double> imagTraceMultiplication(const Lattice<SU3> &L1, const Lattice<SU3> &L2)
{
    /*
     * Function for multiplying lattice and taking the imaginary trace without summing the matrix.
     */
    unsigned long site = 0;
    Lattice<double> latticeSum;
    latticeSum.allocate(L1.m_dim);
    latticeSum.zeros();

    for (unsigned int ix = 0; ix < L1.m_dim[0]; ix++)
    {
        for (unsigned int iy = 0; iy < L1.m_dim[1]; iy++)
        {
            for (unsigned int iz = 0; iz < L1.m_dim[2]; iz++)
            {
                for (unsigned int it = 0; it < L1.m_dim[3]; it++)
                {
                    site = Parallel::Index::getIndex(ix,iy,iz,it);
                    latticeSum[site] += traceImagMultiplication(L1.m_sites[site], L2.m_sites[site]);
                }
            }
        }
    }

    return latticeSum;
}


template <class SU3>
inline double sumRealTrace(const Lattice<SU3> &L)
{
    double latticeSum;
    latticeSum = 0.0;
    for (unsigned long int iSite = 0; iSite < L.m_latticeSize; iSite++)
    {
        latticeSum += L.m_sites[iSite].mat[0] + L.m_sites[iSite].mat[8] + L.m_sites[iSite].mat[16];
    }
    return latticeSum;
}

template <class SU3>
inline double sumRealTraceMultiplication(const Lattice<SU3> &L1, const Lattice<SU3> &L2)
{
    double latticeSum;
    latticeSum = 0.0;
    for (unsigned long int iSite = 0; iSite < L1.m_latticeSize; iSite++)
    {
        latticeSum += traceRealMultiplication(L1.m_sites[iSite],L2.m_sites[iSite]);
    }
    return latticeSum;
}

template <class SU3>
inline Lattice<SU3> inv(Lattice<SU3> &B)
{
    Lattice<SU3> _L;
    _L.allocate(B.m_dim);
    for (unsigned long int iSite = 0; iSite < B.m_latticeSize; iSite++)
    {
        _L.m_sites[iSite] = B.m_sites[iSite].inv();
    }
    return std::move(_L);
}

template <class SU3>
inline Lattice<SU3> inv(Lattice<SU3> &&B)
{
    Lattice<SU3> _L;
    _L.allocate(B.m_dim);
    for (unsigned long int iSite = 0; iSite < B.m_latticeSize; iSite++)
    {
        _L.m_sites[iSite] = B.m_sites[iSite].inv();
    }
    return std::move(_L);
}

template <class SU3>
inline Lattice<SU3> transpose(Lattice<SU3> &B)
{
    Lattice<SU3> _L;
    _L.allocate(B.m_dim);
    for (unsigned long int iSite = 0; iSite < B.m_latticeSize; iSite++)
    {
        _L.m_sites[iSite] = B.m_sites[iSite].transpose();
    }
    return std::move(_L);
}

template <class SU3>
inline Lattice<SU3> transpose(Lattice<SU3> &&B)
{
    Lattice<SU3> _L;
    _L.allocate(B.m_dim);
    for (unsigned long int iSite = 0; iSite < B.m_latticeSize; iSite++)
    {
        _L.m_sites[iSite] = B.m_sites[iSite].transpose();
    }
    return std::move(_L);
}

template <class SU3>
inline Lattice<SU3> conjugate(Lattice<SU3> &B)
{
    Lattice<SU3> _L;
    _L.allocate(B.m_dim);
    for (unsigned long int iSite = 0; iSite < B.m_latticeSize; iSite++)
    {
        _L.m_sites[iSite] = B.m_sites[iSite].conjugate();
    }
    return std::move(_L);
}

template <class SU3>
inline Lattice<SU3> conjugate(Lattice<SU3> &&B)
{
    Lattice<SU3> _L;
    _L.allocate(B.m_dim);
    for (unsigned long int iSite = 0; iSite < B.m_latticeSize; iSite++)
    {
        _L.m_sites[iSite] = B.m_sites[iSite].conjugate();
    }
    return std::move(_L);
}

//////////////////////////////////////////
//////// Communication functions /////////
//////////////////////////////////////////
// Enumerator container for direction
/*!
 * \enum DIR
 *
 * \brief The DIR enum specifies if we are to move BACKWARDS or FORWARDS in the lattice shifts.
 */
enum DIR
{
    BACKWARDS = 0,
    FORWARDS = 1
};

/*!
 * \brief shift copies the lattice L and shifts it in a specified lorentzVector direction DIR.
 *
 * The shift shares a face of the lattice with a neighbouring processsors, and retrieves a corresponding face from a neighboring processor, essentialle shifting the entire lattice one \f$\hat{\mu}\f$.
 *
 * \sa Chapter 5 in https://github.com/hmvege/LQCDMasterThesis contains a detailed description of the sharing process.
 *
 * \param L is the lattiece we are to shift
 * \param direction we are shifting in. Either FORWARDS or BACKWARDS.
 * \param lorentzVector the direction in the Lorentz index \f$\hat{\mu}\f$ we are shifting in. Either 0, 1, 2 or 3.
 * \return a copy of the lattice L that contains the face of a neighboring lattice.
 */
inline Lattice<SU3> shift(const Lattice<SU3> &L, const DIR direction, const unsigned int lorentzVector)
{
    /*
     * Function for shifting lattice in on or another direction.
     * The direction of the matrix is taken care if in what lattice we are passing.
     *
     * The FORWADS/BACKWARDS semantic is quite confusing, but consider this:
     * When shifting e.g. BACKWARDS, we want to update the current lattice at processor P
     * with the lattice from the processor in-front, i.e. shifting it one step backwards,
     *  [x11 x12 x13] [y11 y12 y13]    [x12 x13 y11] [y12 y13 x11]
     *  [x21 x22 x23] [y21 y22 y23] -> [x22 x23 y21] [y22 y23 x21]
     *  [x31 x32 x33] [y31 y32 y33]    [x32 x33 y31] [y32 y33 x31]
     * The shifted lattice, is then returned.
     */
    Lattice<SU3> _L;
    _L.allocate(L.m_dim);// MOVE THIS TO INITIALIZATION/HEADER-THING?
    std::vector<SU3> sendCube; // Move indexes to index in order to avoid 2 integer multiplications)
    std::vector<SU3> recvCube; // MOVE THIS TO HEADER; SO WE DONT ALLOCATE EVERY TIME!
    // INSTEAD OF CUBE INDEX, JUST DO INDEX WITH DIMENSION SET TO ZERO
    MPI_Request sendReq,recvReq;
    switch(direction) {
    case BACKWARDS: {
        switch(lorentzVector) {
        case 0: {
            sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]); // Four of these can actually be stored globally
            recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            /* Max memory usage: 48*48*48*96 /w 512 procs -->  12 12 12 12 --> 4 cubes of size 12^3 = 1728*18 bytes
             * --> 1728*28 / 1024(to kilobytes) / 1024(to megabytes) = 0.03 MB per processor.
             * Maximum use of 4 volumes, one for each direction(assuming that spatial directionality may vary) --> 0.12 MB in total for this part
             */
            // Populates package to send
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) { // Ensure cube indexes match before and after in order to map correctly!!
                        sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(L.m_dim[0]-1,iy,iz,it)];
                    }
                }
            }
            // Sends and receives packages
            MPI_Isend(&sendCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(1),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(0),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 1; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                           _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix-1,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received, then populates face cube with remaining results.
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            // Repopulates the lattice with missing cube
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                       _L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)] = recvCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 1: {
            sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,L.m_dim[1]-1,iz,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(3),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(2),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 1; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy-1,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)] = recvCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 2: {
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,L.m_dim[2]-1,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(5),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(4),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 1; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz-1,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,0,it)] = recvCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        case 3: {
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,L.m_dim[3]-1)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(7),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(6),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 1; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it-1)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                       _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,0)] = recvCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        }
        break;
    }
    case FORWARDS: {
        switch(lorentzVector) {
        case 0: { // x direction
            sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)];
                    }
                }
            }
            // Sends and receives packages
            MPI_Isend(&sendCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(0),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(1),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0] - 1; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                           _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix+1,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received, then populates face cube with remaining results.
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                       _L.m_sites[Parallel::Index::getIndex(L.m_dim[0]-1,iy,iz,it)] = recvCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 1: { // y direction
            sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(2),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(3),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1] - 1; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                           _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy+1,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                       _L.m_sites[Parallel::Index::getIndex(ix,L.m_dim[1]-1,iz,it)] = recvCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 2: { // z direction
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,0,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(4),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(5),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2] - 1; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz+1,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,L.m_dim[2]-1,it)] = recvCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        case 3: { // t direction
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,0)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(6),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(7),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3] - 1; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it+1)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,L.m_dim[3]-1)] = recvCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        }
        break;
    }
    }
    return _L;
}

inline Lattice<SU3> shift(Lattice<SU3> &&L, const DIR direction, const unsigned int lorentzVector)
{
    /*
     * Function for shifting lattice in on or another direction.
     * The direction of the matrix is taken care if in what lattice we are passing.
     *
     * The FORWADS/BACKWARDS semantic is quite confusing, but consider this:
     * When shifting e.g. BACKWARDS, we want to update the current lattice at processor P
     * with the lattice from the processor in-front, i.e. shifting it one step backwards,
     *  [x11 x12 x13] [y11 y12 y13]    [x12 x13 y11] [y12 y13 x11]
     *  [x21 x22 x23] [y21 y22 y23] -> [x22 x23 y21] [y22 y23 x21]
     *  [x31 x32 x33] [y31 y32 y33]    [x32 x33 y31] [y32 y33 x31]
     * The shifted lattice, is then returned.
     */
    Lattice<SU3> _L;
    _L.allocate(L.m_dim);// MOVE THIS TO INITIALIZATION/HEADER-THING?
    std::vector<SU3> sendCube; // Move indexes to index in order to avoid 2 integer multiplications)
    std::vector<SU3> recvCube; // MOVE THIS TO HEADER; SO WE DONT ALLOCATE EVERY TIME!
    // INSTEAD OF CUBE INDEX, JUST DO INDEX WITH DIMENSION SET TO ZERO
    MPI_Request sendReq,recvReq;
    switch(direction) {
    case BACKWARDS: {
        switch(lorentzVector) {
        case 0: {
            sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]); // Four of these can actually be stored globally
            recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            /* Max memory usage: 48*48*48*96 /w 512 procs -->  12 12 12 12 --> 4 cubes of size 12^3 = 1728*18 bytes
             * --> 1728*28 / 1024(to kilobytes) / 1024(to megabytes) = 0.03 MB per processor.
             * Maximum use of 4 volumes, one for each direction(assuming that spatial directionality may vary) --> 0.12 MB in total for this part
             */
            // Populates package to send
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) { // Ensure cube indexes match before and after in order to map correctly!!
                        sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(L.m_dim[0]-1,iy,iz,it)];
                    }
                }
            }
            // Sends and receives packages
            MPI_Isend(&sendCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(1),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(0),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 1; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                           _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix-1,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received, then populates face cube with remaining results.
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            // Repopulates the lattice with missing cube
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                       _L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)] = recvCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 1: {
            sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,L.m_dim[1]-1,iz,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(3),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(2),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 1; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy-1,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)] = recvCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 2: {
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,L.m_dim[2]-1,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(5),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(4),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 1; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz-1,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,0,it)] = recvCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        case 3: {
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,L.m_dim[3]-1)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(7),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(6),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 1; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it-1)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                       _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,0)] = recvCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        }
        break;
    }
    case FORWARDS: {
        switch(lorentzVector) {
        case 0: { // x direction
            sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)];
                    }
                }
            }
            // Sends and receives packages
            MPI_Isend(&sendCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(0),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(1),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0] - 1; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                           _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix+1,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received, then populates face cube with remaining results.
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                       _L.m_sites[Parallel::Index::getIndex(L.m_dim[0]-1,iy,iz,it)] = recvCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 1: { // y direction
            sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(2),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(3),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1] - 1; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                           _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy+1,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                       _L.m_sites[Parallel::Index::getIndex(ix,L.m_dim[1]-1,iz,it)] = recvCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])];
                    }
                }
            }
            return _L;
        }
        case 2: { // z direction
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,0,it)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(4),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Neighbours::get(5),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2] - 1; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz+1,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,L.m_dim[2]-1,it)] = recvCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        case 3: { // t direction
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            // Populates package to send
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,0)];
                    }
                }
            }
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(6),0,Parallel::ParallelParameters::ACTIVE_COMM,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Neighbours::get(7),0,Parallel::ParallelParameters::ACTIVE_COMM,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (unsigned int it = 0; it < L.m_dim[3] - 1; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it+1)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (unsigned int ix = 0; ix < L.m_dim[0]; ix++) {
                for (unsigned int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (unsigned int iz = 0; iz < L.m_dim[2]; iz++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,L.m_dim[3]-1)] = recvCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            return _L;
        }
        }
        break;
    }
    }
    return _L;
}

#endif // LATTICE_H
