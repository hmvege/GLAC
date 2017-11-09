#ifndef NEIGHBOURLIST_H
#define NEIGHBOURLIST_H

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

struct NeighbourList
{
    NeighbourList();
    ~NeighbourList();
    // Neighbour list:
    int list[8];

    int rank;
    void print();
};

#endif // NEIGHBOURLIST_H
