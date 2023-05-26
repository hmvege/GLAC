#include "neighbourlist.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

NeighbourList:: NeighbourList()
{
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */
}

NeighbourList::~NeighbourList()
{
}

/*!
 * \brief NeighbourList::print prints the geometry of a single processor.
 *
 * Usefull for debugging.
 */
void NeighbourList::print()
{
    cout << "\nProcess rank = " << rank << endl;
    cout << "x-1 = " << list[0] << " x+1 = " << list[1] << endl;
    cout << "y-1 = " << list[2] << " y+1 = " << list[3] << endl;
    cout << "z-1 = " << list[4] << " z+1 = " << list[5] << endl;
    cout << "t-1 = " << list[6] << " t+1 = " << list[7] << "\n" << endl;
}
