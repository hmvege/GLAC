/*!
 * \brief
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef NEIGHBOURLIST_H
#define NEIGHBOURLIST_H

struct NeighbourList
{
    NeighbourList();
    ~NeighbourList();
    // Neighbour list:
    int list[8];

    int rank;
    void print();

    int &operator [](int i) { return list[i]; }
};

#endif // NEIGHBOURLIST_H
