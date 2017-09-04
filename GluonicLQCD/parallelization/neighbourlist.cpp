#include "neighbourlist.h"

// TEMP DIAGNOSTICS
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
    list = new int[8];
}

NeighbourList::~NeighbourList()
{
    cout << "DELETING" << endl;
    delete [] list;
}

void NeighbourList::print()
{
//    for (int i = 0; i < 8; i++) {
////        if (list.at(i) != NULL) {
////            cout << "Segfault at i = " << i << endl;
//            cout << "list["<<i<<"] = "<<list[i]<<endl;
////        }
//    }
    cout << "\nProcess rank = " << rank << endl;
    cout << "x-1 = " << list[0] << " x+1 = " << list[1] << endl;
    cout << "y-1 = " << list[2] << " y+1 = " << list[3] << endl;
    cout << "z-1 = " << list[4] << " z+1 = " << list[5] << endl;
    cout << "t-1 = " << list[6] << " t+1 = " << list[7] << "\n" << endl;
}
