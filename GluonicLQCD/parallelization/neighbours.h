#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include <iostream>

#include "neighbourlist.h"

typedef int (*indexCubeFinderArray) (int n1, int n2, int n3);

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

    // Contigious index finder for cubes
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */
//    static int indexCube0(int n1, int n2, int n3) { return (m_Nt*(m_Nz*n1 + n2) + n3); } // x-1, y,z,t
//    int indexCube1(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*m_Nx + n1) + n2) + n3); } // x+1, y,z,t
//    int indexCube2(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*n1) + n2) + n3); } // y-1, x,z,t
//    int indexCube3(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*n1 + m_Ny) + n2) + n3); } // y+1, x,z,t
//    int indexCube4(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*n1 + n2)) + n3); } // z-1, x,y,t
//    int indexCube5(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*n1 + n2) + m_Nz) + n3); } // z+1, x,y,t
//    int indexCube6(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*n1 + n2) + n3)); } // t-1, x,y,z
//    int indexCube7(int n1, int n2, int n3) { return (m_Nt*(m_Nz*(m_Ny*n1 + n2) + n3) + m_Nz); } // t+1, x,y,z
//    void generateCubeIndexFunctionList();
//    void generateCubeIndexes();

    // Functions for finding correct direction of neighbouring lattices
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
        if ((((Np / m_Nx) + 1) % m_Ny) == 0) {
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
            int x_count = Np % m_Nx;
            int y_count = ((Np/m_Nx) % m_Ny) * m_Nx;
            int t_count = ((Np/(m_Nx*m_Ny*m_Nz)) % m_Nz) * (m_Nx*m_Ny*m_Nt);
//            return Np % m_Nx + ((Np/m_Nx) % m_Ny) * m_Nx;
            return x_count + y_count + t_count;
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
    Neighbours();
    ~Neighbours();

    void initialize(int processRank, int numproc, int * processorsPerDim);

//    // Index function for specific face
//    indexCubeFinderArray *cubeIndexFunctions;
//    int ** cubeIndex;

    NeighbourList* getNeighbours(int Np);
    int getListLength() { return m_numproc; }
};

#endif // NEIGHBOURS_H
