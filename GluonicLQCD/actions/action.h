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
    int *m_N;
    // For handling the shift-method in parallelization
    std::vector<int> indexes;
    IndexOrganiser *m_Index = nullptr;
public:
    Action();
//    Action(int *N);
    virtual ~Action();
    virtual double getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu);
    virtual void computeStaple(Links *lattice, int i, int j, int k, int l, int mu);
    // Setters
    void setN(int *N);
    void initializeIndexHandler(IndexOrganiser *Index);
};

#endif // ACTION_H
