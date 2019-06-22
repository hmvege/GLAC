/*!
 * \class Neighbours
 *
 * \brief Class for setting up the NeighbourList and stores them.
 *
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "neighbourlist.h"
#include <vector>

namespace Parallel {
class Neighbours
{
private:
    static int m_processRank;
    static int m_numproc;
    static int m_Nx, m_Ny, m_Nz, m_Nt; // Prosessors per dimension
    static long m_P[4]; // Prosessor coordinate
    static void generateNeighbourList();
    static std::vector<NeighbourList> m_neighbourLists;

    // Contigious index finder for cubes
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */

    // Functions for finding correct direction of neighbouring lattices
    static inline int getXPlusOne(const int Np)
    {
        if (((Np + 1) % m_Nx) == 0)
        {
            return (Np/m_Nx) * m_Nx;
        }
        else
        {
            return Np + 1;
        }
    }

    static inline int getXMinusOne(const int Np)
    {
        int x_count = (Np - 1 + m_Nx) % m_Nx;
        int y_count = (Np/m_Nx) * m_Nx;
        return x_count + y_count;
    }

    static inline int getYPlusOne(const int Np)
    {
        if ((((Np / m_Nx) + 1) % m_Ny) == 0)
        {
            int x_count = Np % m_Nx;
            int y_count = (Np / (m_Ny*m_Nx)) * (m_Nx*m_Ny);
            return x_count + y_count;
        }
        else
        {
            return Np + m_Nx;
        }
    }

    static inline int getYMinusOne(const int Np)
    {
        int x_count = Np % m_Nx;
        int y_count = (((Np/m_Nx) - 1 + m_Ny) % m_Ny) * m_Nx;
        int z_count = (Np/(m_Nx * m_Ny)) * (m_Nx*m_Ny);
        return x_count + y_count + z_count;
    }

    static inline int getZPlusOne(const int Np)
    {
        if (((Np/(m_Nx*m_Ny) + 1) % m_Nz) == 0)
        {
            int x_count = Np % m_Nx;
            int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
            int t_count = ((Np/(m_Nx*m_Ny*m_Nz)) % m_Nt) * (m_Nx*m_Ny*m_Nz);
            return x_count + y_count + t_count;
        }
        else
        {
            return Np + m_Nx*m_Ny;
        }
    }

    static inline int getZMinusOne(const int Np)
    {
        int x_count = Np % m_Nx;
        int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
        int z_count = ((Np/(m_Nx*m_Ny) - 1 + m_Nz) % m_Nz) * (m_Nx*m_Ny);
        int t_count = (Np/(m_Nx*m_Ny*m_Nz)) * (m_Nx*m_Ny*m_Nz);
        return x_count + y_count + z_count + t_count;
    }

    static inline int getTPlusOne(const int Np)
    {
        if ((((Np/(m_Nx*m_Ny*m_Nz)) + 1) % m_Nt) == 0)
        {
            int x_count = Np % m_Nx;
            int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
            int z_count = (Np/(m_Nx*m_Ny) % m_Nz) * (m_Nx*m_Ny);
            return x_count + y_count + z_count;
        }
        else
        {
            return Np + m_Nx*m_Ny*m_Nz;
        }
    }

    static inline int getTMinusOne(const int Np)
    {
        int x_count = Np % m_Nx;
        int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
        int z_count = (Np/(m_Nx*m_Ny) % m_Nz) * (m_Nx*m_Ny);
        int t_count = ((Np/(m_Nx*m_Ny*m_Nz) - 1 + m_Nt) % m_Nt) * (m_Nx*m_Ny*m_Nz);
        return x_count + y_count + z_count + t_count;
    }
public:
    Neighbours();
    ~Neighbours();

    static void initialize(const int processRank, const int numproc, int * processorsPerDim);

    /*!
     * \brief get returns the neighbourlist for the the calling processor/rank and the direction given py iProcDir.
     * \param iProcDir the direction in which to look up its neighbour.
     * \return the rank of its neighbour.
     */
    static int get(const int iProcDir) { return m_neighbourLists[m_processRank][iProcDir]; } // Returns neighbour list for iP processor.

    // Getters
    static NeighbourList* getNeighbours(const int Np);
    static int getListLength() { return m_numproc; }
    static long getProcessorDimensionPosition(const int dim) { return m_P[dim]; }
};
}

#endif // NEIGHBOURS_H
