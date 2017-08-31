#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include <iostream>

#include "neighbourlist.h"

class Neighbours
{
private:
    int m_processRank;
    int m_numproc;
    int m_Nx;
    int m_Ny;
    int m_Nz;
    int m_Nt;
    void generateNeighbourList();
    NeighbourList * neighbourLists;

    inline int getXPlusOne(int Np) {
        if (((Np + 1) % m_Nx) == 0) {
            return (Np/m_Nx) * m_Nx;
        } else {
            return Np + 1;
        }
    }

    inline int getXMinusOne(int Np) {
        int x_count = (Np - 1 + m_Nx) % m_Nx;
        int y_count = (Np/m_Nx) * m_Nx;
        return x_count + y_count;
    }

    inline int getYPlusOne(int Np) {
        if ((((Np / m_Nx) + 1) % m_Nx) == 0) {
            int x_count = Np % m_Nx;
            int y_count = (Np / (m_Ny*m_Nx)) * (m_Nx*m_Ny);
            return x_count + y_count;
        } else {
            return Np + m_Nx;
        }
    }

    inline int getYMinusOne(int Np) {
        int x_count = Np % m_Nx;
        int y_count = (((Np/m_Nx) - 1 + m_Ny) % m_Ny) * m_Nx;
        int z_count = (Np/(m_Nx * m_Ny)) * (m_Nx*m_Ny);
        return x_count + y_count + z_count;
    }

    inline int getZPlusOne(int Np) {
        if (((Np/(m_Nx*m_Ny) + 1) % m_Nz) == 0) {
            return Np % m_Nx + ((Np/m_Nx) % m_Ny) * m_Nx;
        } else {
            return Np + m_Nx*m_Ny;
        }
    }

    inline int getZMinusOne(int Np) {
        int x_count = Np % m_Nx;
        int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
        int z_count = ((Np/(m_Nx*m_Ny) - 1 + m_Nz) % m_Nz) * (m_Nx*m_Ny);
        int t_count = (Np/(m_Nx*m_Ny*m_Nz)) * (m_Nx*m_Ny*m_Nz);
        return x_count + y_count + z_count + t_count;
    }

    inline int getTPlusOne(int Np) {
        if ((((Np/(m_Nx*m_Ny*m_Nz)) + 1) % m_Nt) == 0) {
            int x_count = Np % m_Nx;
            int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
            int z_count = (Np/(m_Nx*m_Ny) % m_Nz) * (m_Nx*m_Ny);
            return x_count + y_count + z_count;
        } else {
            return Np + m_Nx*m_Ny*m_Nz;
        }
    }

    inline int getTMinusOne(int Np) {
        int x_count = Np % m_Nx;
        int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
        int z_count = (Np/(m_Nx*m_Ny) % m_Nz) * (m_Nx*m_Ny);
        int t_count = ((Np/(m_Nx*m_Ny*m_Nz) - 1 + m_Nt) % m_Nt) * (m_Nx*m_Ny*m_Nz);
        return x_count + y_count + z_count + t_count;
    }
public:
    Neighbours(int processRank, int numproc, int *processorsPerDim);
    ~Neighbours();

    NeighbourList getNeighbours(int Np) { return neighbourLists[Np]; }
};

#endif // NEIGHBOURS_H
