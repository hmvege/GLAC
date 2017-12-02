#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <mpi.h>
#include "math/complex.h"
#include "parallelization/index.h"
#include "parallelization/communicator.h"

template <class T>
class Lattice
{
public:
    std::vector<T> m_sites;
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
    // Lattice value setters
    void identity();
    void zeros();
    Lattice<T> inv();
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
///////// Lattice value setters //////////
//////////////////////////////////////////
// Value setters of the lattice
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

template <class T>
inline Lattice<T> Lattice<T>::inv() {
    for (int iSite = 0; iSite < m_latticeSize; iSite++) {
        m_sites[iSite].inv();
    }
    return *this; // Okay to do this?? Do I need to initialize a new lattice??
}

// Allocates memory to the lattice. Has to be called every time(unless we are copying)
template <class T>
inline void Lattice<T>::allocate(std::vector<int> dim) {
    m_latticeSize = dim[0] * dim[1] * dim[2] * dim[3];
    m_dim = dim;
    m_sites.resize(m_latticeSize);
}

//////////////////////////////////////////
////////// Lattice functions /////////////
//////////////////////////////////////////
template <class T>
inline Lattice<double> realTrace(Lattice<T> L)
{
    Lattice<double> tempTraceSum(L.m_dim);
    for (int iSite = 0; iSite < L.m_latticeSize; iSite++) {
        tempTraceSum[iSite] = (L[iSite][0] + L[iSite][8] + L[iSite][16]);
    }
    return tempTraceSum;
}

template <class T>
inline Lattice<double> imagTrace(Lattice<T> L)
{
    Lattice<double> tempTraceSum(L.m_dim);
    for (int iSite = 0; iSite < L.m_latticeSize; iSite++) {
        tempTraceSum[iSite] = (L[iSite][1] + L[iSite][9] + L[iSite][17]);
    }
    return tempTraceSum;
}

template <class T>
inline Lattice<complex> trace(Lattice<T> L)
{
    Lattice<complex> tempTraceSum(L.m_dim);
    for (int iSite = 0; iSite < L.m_latticeSize; iSite++) {
        tempTraceSum[iSite] = complex(L[iSite][0] + L[iSite][8] + L[iSite][16],
                L[iSite][1] + L[iSite][9] + L[iSite][17]);
    }
    return tempTraceSum;
}

template <class T>
inline T sum(Lattice<T> L)
{
    T latticeSum;
    latticeSum = 0.0;
    for (int iSite = 0; iSite < L.m_latticeSize; iSite++) {
        latticeSum += L[iSite];
    }
    return latticeSum;
}

//////////////////////////////////////////
//////// Communication functions /////////
//////////////////////////////////////////
// Enumerator container for direction
enum DIR {
    BACKWARDS = 0,
    FORWARDS = 1
};

//template <class T>
inline Lattice<SU3> shift(Lattice<SU3> L, DIR direction, int lorentzVector)
{
    /*
     * Function for shifting lattice in on or another direction.
     * The direction of the matrix is taken care if in what lattice we are passing.
     */
    Lattice<SU3> _L;
    _L.allocate(L.m_dim);// MOVE THIS TO INITIALIZATION/HEADER-THING?
    std::vector<SU3> sendCube; // Move indexes to index in order to avoid 2 integer multiplications)
    std::vector<SU3> recvCube; // MOVE THIS TO HEADER; SO WE DONT ALLOCATE EVERY TIME!
    MPI_Request sendReq,recvReq;

    switch(direction) {
    case BACKWARDS: {
        switch(lorentzVector) {
        case 0: {
            sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (int iy = 0; iy < L.m_dim[1]; iy++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) { // Ensure cube indexes match before and after in order to map correctly!!
                        sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(L.m_dim[0]-1,iy,iz,it)];
                    }
                }
            }
            // Sends and receives packages
            MPI_Isend(&sendCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[1],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[0],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 1; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received, then populates face cube with remaining results.
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            // Repopulates the lattice with missing cube
            for (int iy = 0; iy < L.m_dim[1]; iy++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)] = recvCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                    }
                }
            }
            break;
        }
        case 1: {
            sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,L.m_dim[1]-1,iz,it)];
                    }
                }
            }
            // MPI_Request req;
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[3],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[2],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 1; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)] = recvCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])];
                    }
                }
            }
            break;
        }
        case 2: {
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            // Populates package to send
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,L.m_dim[2]-1,it)];
                    }
                }
            }
            // MPI_Request req;
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[5],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[4],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 1; iz < L.m_dim[2]; iz++) {
                        for (int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,0,it)] = recvCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            break;
        }
        case 3: {
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            // Populates package to send
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,L.m_dim[3]-1)];
                    }
                }
            }
            // MPI_Request req;
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Communicator::m_NLists[7],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Communicator::m_NLists[6],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (int it = 1; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,0)] = recvCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            break;
        }
        }
        break;
    }
    case FORWARDS: {
        switch(lorentzVector) {
        case 0: { // x direction
            sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]); // Four of these can actually be stored globally
            recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);
            /* Max memory usage: 48*48*48*96 /w 512 procs -->  12 12 12 12 --> 4 cubes of size 12^3 = 1728*18 bytes
             * --> 1728*28 / 1024(to kilobytes) / 1024(to megabytes) = 0.03 MB.
             * Maximum use of 4 volumes, one for each direction(assuming that spatial directionality may vary) --> 0.12 MB in total for this part
             */
            // Populates package to send
            for (int iy = 0; iy < L.m_dim[1]; iy++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) { // Ensure cube indexes match before and after in order to map correctly!!
                        sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)];
                    }
                }
            }
            // Sends and receives packages
            MPI_Isend(&sendCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[0],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[1],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0] - 1; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received, then populates face cube with remaining results.
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            for (int iy = 0; iy < L.m_dim[1]; iy++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(L.m_dim[0]-1,iy,iz,it)] = recvCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                    }
                }
            }
            break;
        }
        case 1: { // y direction
            sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
            // Populates package to send
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)];
                    }
                }
            }
            // MPI_Request req;
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[2],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[3],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1] - 1; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,L.m_dim[1]-1,iz,it)] = recvCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])];
                    }
                }
            }
            break;
        }
        case 2: { // z direction
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]);
            // Populates package to send
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,0,it)];
                    }
                }
            }
            // MPI_Request req;
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[4],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[5],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2] - 1; iz++) {
                        for (int it = 0; it < L.m_dim[3]; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,L.m_dim[2]-1,it)] = recvCube[Parallel::Index::cubeIndex(ix,iy,it,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            break;
        }
        case 3: { // t direction
            sendCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            recvCube.resize(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]);
            // Populates package to send
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        sendCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,0)];
                    }
                }
            }
            // MPI_Request req;
            MPI_Isend(&sendCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Communicator::m_NLists[6],0,MPI_COMM_WORLD,&sendReq);
            MPI_Irecv(&recvCube.front(),18*(L.m_dim[0]*L.m_dim[1]*L.m_dim[2]),MPI_DOUBLE,Parallel::Communicator::m_NLists[7],0,MPI_COMM_WORLD,&recvReq);
            // Populates shifted lattice by elements not required to share
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        for (int it = 0; it < L.m_dim[3] - 1; it++) {
                            _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                        }
                    }
                }
            }
            // Ensures all results have been sent and received
            MPI_Wait(&recvReq,MPI_STATUS_IGNORE);
            MPI_Wait(&sendReq,MPI_STATUS_IGNORE);
            for (int ix = 0; ix < L.m_dim[0]; ix++) {
                for (int iy = 0; iy < L.m_dim[1]; iy++) {
                    for (int iz = 0; iz < L.m_dim[2]; iz++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,L.m_dim[3]-1)] = recvCube[Parallel::Index::cubeIndex(ix,iy,iz,L.m_dim[0],L.m_dim[1])];
                    }
                }
            }
            break;
        }
        }
        break;
    }
    }
    return _L; // Should never need to go
}

template <class T>
inline Lattice<T> _lorentzSwitch(Lattice<T> L, DIR direction, int lorentzVector)
{
    Lattice<T> _L;
    std::vector<T> sendCube; // Move indexes to index in order to avoid 2 integer multiplications)
    std::vector<T> recvCube; // MOVE THIS TO HEADER; SO WE DONT ALLOCATE EVERY TIME!
    MPI_Request req;
    switch(lorentzVector) {
    case 0: // x direction
        sendCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]); // Four of these can actually be stored globally
        /* Max memory usage: 48*48*48*96 /w 512 procs -->  12 12 12 12 --> 4 cubes of size 12^3 = 1728*18 bytes
         * --> 1728*28 / 1024(to kilobytes) / 1024(to megabytes) = 0.03 MB.
         * Maximum use of 4 volumes, one for each direction(assuming that spatial directionality may vary) --> 0.12 MB in total for this part
         */
        recvCube.resize(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]);

        for (int iy = 0; iy < L.m_dim[1]; iy++) {
            for (int iz = 0; iz < L.m_dim[2]; iz++) {
                for (int it = 0; it < L.m_dim[3]; it++) {
                    sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(0,iy,iz,it)];
                }
            }
        }
        MPI_Isend(&sendCube,18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[direction],0,MPI_COMM_WORLD,&req);
        MPI_Irecv(&recvCube,18*(L.m_dim[1]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[abs(direction - 1)],0,MPI_COMM_WORLD,&req);

        for (int ix = 1; ix < L.m_dim[0] - 1; ix++) {
            for (int iy = 0; iy < L.m_dim[1]; iy++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                    }
                }
            }
        }
        MPI_Wait(&req,MPI_STATUS_IGNORE);
        for (int ix = 1; ix < L.m_dim[0] - 1; ix++) {
            for (int iy = 0; iy < L.m_dim[1]; iy++) {
                for (int iz = 0; iz < L.m_dim[2]; iz++) {
                    for (int it = 0; it < L.m_dim[3]; it++) {
                        _L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)] = L.m_sites[Parallel::Index::getIndex(ix,iy,iz,it)];
                    }
                }
            }
        }
        for (int iy = 0; iy < L.m_dim[1]; iy++) {
            for (int iz = 0; iz < L.m_dim[2]; iz++) {
                for (int it = 0; it < L.m_dim[3]; it++) {
                    _L.m_sites[Parallel::Index::getIndex(L.m_dim[1]-1,iy,iz,it)] = sendCube[Parallel::Index::cubeIndex(iy,iz,it,L.m_dim[1],L.m_dim[2])];
                }
            }
        }

    case 1: // y direction
        //        std::vector<T> sendCube; // Move indexes to index in order to avoid 2 integer multiplications)
        //        std::vector<T> recvCube; // MOVE THIS TO HEADER; SO WE DONT ALLOCATE EVERY TIME!
        //        sendCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);
        //        recvCube.resize(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]);

        //        for (int ix = 0; ix < L.m_dim[0]; ix++) {
        //            for (int iz = 0; iz < L.m_dim[2]; iz++) {
        //                for (int it = 0; it < L.m_dim[3]; it++) {
        //                    sendCube[Parallel::Index::cubeIndex(ix,iz,it,L.m_dim[0],L.m_dim[2])] = L.m_sites[Parallel::Index::getIndex(ix,0,iz,it)];
        //                }
        //            }
        //        }
        ////        MPI_Request req;
        //        MPI_Isend(&sendCube,18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[direction],0,MPI_COMM_WORLD,&req);
        //        MPI_Irecv(&recvCube,18*(L.m_dim[0]*L.m_dim[2]*L.m_dim[3]),MPI_DOUBLE,Parallel::Communicator::m_NLists[abs(direction - 1)],0,MPI_COMM_WORLD,&req);
    case 2: // z direction
        ;
    case 3: // t direction
        ;
    }
    return _L;
}

#endif // LATTICE_H
