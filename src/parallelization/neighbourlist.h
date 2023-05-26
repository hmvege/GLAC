/*!
 * \class NeighbourLists
 *
 * \brief Class for storing the nearest closest neighbours of a processors.
 *
 * Stores a static int array of 8 elements, since there is a neighbour in each x, y, z and t direction.
 *
 * Neighbour list values defined as:
 * 0: x-1 | 1: x+1
 * 2: y-1 | 3: y+1
 * 4: z-1 | 5: z+1
 * 6: t-1 | 7: t+1
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

    //! Neighbour list
    int list[8];

    //! Rank of processor holding NeighbourList.
    int rank;
    void print();

    /*!
     * \brief operator [] overloaded list accessing tool.
     * \param i direction to look up.
     * \return the process rank in given i direction.
     */
    int &operator [](int i) { return list[i]; }
};

#endif // NEIGHBOURLIST_H
