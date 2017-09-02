#include "neighbours.h"
#include "neighbourlist.h"

// DEBUGGING
#include <iostream>

//Neighbours::Neighbours(int processRank, int numproc, int * processorsPerDim)
Neighbours::Neighbours()
{
//    m_processRank = processRank;
//    m_numproc = numproc;
//    m_Nx = processorsPerDim[0];
//    m_Ny = processorsPerDim[1];
//    m_Nz = processorsPerDim[2];
//    m_Nt = processorsPerDim[3];
//    neighbourLists = new NeighbourList[m_numproc];
//    generateNeighbourList();
}

Neighbours::~Neighbours()
{
    delete [] neighbourLists;
//    for (int i = 0; i < 4; i++) {
//        delete [] cubeIndex[i];
//    }
//    delete [] cubeIndex;
//    delete [] cubeIndexFunctions;
}

void Neighbours::initialize(int processRank, int numproc, int * processorsPerDim) {
    m_processRank = processRank;
    m_numproc = numproc;
    m_Nx = processorsPerDim[0];
    m_Ny = processorsPerDim[1];
    m_Nz = processorsPerDim[2];
    m_Nt = processorsPerDim[3];
    neighbourLists = new NeighbourList[m_numproc];
    generateNeighbourList();
//    generateCubeIndexFunctionList();
//    generateCubeIndexes();
}

//void Neighbours::generateCubeIndexes()
//{
//    cubeIndex = new int*[4];
//    for (int i = 0; i < 4; i++) {
//        cubeIndex[i] = new int[3];
//    }
//    // x-direction
//    cubeIndex[0][0] = m_Ny;
//    cubeIndex[0][1] = m_Nz;
//    cubeIndex[0][2] = m_Nt;
//    // y-direction
//    cubeIndex[1][0] = m_Nx;
//    cubeIndex[1][1] = m_Nz;
//    cubeIndex[1][2] = m_Nt;
//    // y-direction
//    cubeIndex[2][0] = m_Nx;
//    cubeIndex[2][1] = m_Ny;
//    cubeIndex[2][2] = m_Nt;
//    // y-direction
//    cubeIndex[3][0] = m_Nx;
//    cubeIndex[3][1] = m_Ny;
//    cubeIndex[3][2] = m_Nz;
//}

//void Neighbours::generateCubeIndexFunctionList()
//{
//    cubeIndexFunctions = new indexCubeFinderArray[8];
//    cubeIndexFunctions[0] = &indexCube0;
//    cubeIndexFunctions[1] = &indexCube1;
//    cubeIndexFunctions[2] = &indexCube2;
//    cubeIndexFunctions[3] = &indexCube3;
//    cubeIndexFunctions[4] = &indexCube4;
//    cubeIndexFunctions[5] = &indexCube5;
//    cubeIndexFunctions[6] = &indexCube6;
//    cubeIndexFunctions[7] = &indexCube7;
//}

void Neighbours::generateNeighbourList()
{
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */
    for (int Np = 0; Np < m_numproc; Np++) {
        neighbourLists[Np].rank = Np;
        neighbourLists[Np].list[0] = getXMinusOne(Np);
        neighbourLists[Np].list[1] = getXPlusOne(Np);
        neighbourLists[Np].list[2] = getYMinusOne(Np);
        neighbourLists[Np].list[3] = getYPlusOne(Np);
        neighbourLists[Np].list[4] = getZMinusOne(Np);
        neighbourLists[Np].list[5] = getZPlusOne(Np);
        neighbourLists[Np].list[6] = getTMinusOne(Np);
        neighbourLists[Np].list[7] = getTPlusOne(Np);
    }
}

NeighbourList* Neighbours::getNeighbours(int Np) {
    /*
     * Must return a list. Therefore, returning a reference. Failure to do so, will result in segfault error.
     * Arguments:
     *  Np  : process rank
     */
    return &neighbourLists[Np];
}
