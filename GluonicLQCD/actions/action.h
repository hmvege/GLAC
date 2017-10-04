#ifndef ACTION_H
#define ACTION_H

#include "links.h"
#include "functions.h"
#include "parallelization/indexorganiser.h"
#include "parallelization/neighbours.h"
#include <vector>

class Action
{
protected:
    // Index length array
    unsigned int *m_N;
    // For handling the shift-method in parallelization
    std::vector<int> indexes;
    IndexOrganiser *m_Index = nullptr;
public:
    Action();
//    Action(int *N);
    virtual ~Action();
    virtual double getDeltaAction(Links * lattice, SU3 U, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    virtual void computeStaple(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    // Setters
    void setN(unsigned int *N);
    void initializeIndexHandler(IndexOrganiser *Index);

//    inline void lorentzIndex(int mu, int *lorentzIndices)
//    {
//        /*
//         * Fills the lorentz array with correct indices
//         */
//        for (int i = 0; i < 4; i++)
//        {
//            lorentzIndices[i] = 0;
//        }
//        lorentzIndices[mu] = 1;
//    }
};



#endif // ACTION_H
