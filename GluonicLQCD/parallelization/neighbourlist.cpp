#include "neighbourlist.h"

// TEMP DIAGNOSTICS
#include <iostream>
#include <iomanip>

NeighbourList::NeighbourList()
{
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */
    list = new int[8];
}

NeighbourList::~NeighbourList()
{
    delete [] list;
}

void NeighbourList::print()
{
    std::cout << "\nProcess rank = " << rank << std::endl;
    std::cout << "x-1 = " << list[0] << " x+1 = " << list[1] << std::endl;
    std::cout << "y-1 = " << list[2] << " y+1 = " << list[3] << std::endl;
    std::cout << "z-1 = " << list[4] << " z+1 = " << list[5] << std::endl;
    std::cout << "t-1 = " << list[6] << " t+1 = " << list[7] << "\n" << std::endl;
}
