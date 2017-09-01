#ifndef NEIGHBOURLIST_H
#define NEIGHBOURLIST_H


class NeighbourList
{
public:
    NeighbourList();
    ~NeighbourList();
    // Neighbour list:
    int *list;

//    int *dimensionsToCubesToShareWithNeighbour;
    int rank;
    void print();
};

#endif // NEIGHBOURLIST_H
