#include "neighbourlist.h"

NeighbourList::NeighbourList()
{
    /*
     * Neighbour list values defined as:
     * 0: x-1 | 1: x+1
     * 2: y-1 | 3: y+1
     * 4: z-1 | 5: z+1
     * 6: t-1 | 7: t+1
     */
    m_list = new int[8];
}

NeighbourList::~NeighbourList()
{
    delete [] m_list;
}
