#ifndef ACTION_H
#define ACTION_H

#include "links.h"
#include "functions.h"
#include "parallelization/indexorganiser.h"
#include <vector>

class Action
{
protected:
    // Index length array
    int *m_N;
    // Lorentz indices arrays
    int *muIndex;
    int *nuIndex;
    // For handling the shift-method in parallelization
    std::vector<int> indexes;
    IndexOrganiser *m_Index;
public:
    Action();
    Action(int *N);
    virtual ~Action();
    virtual double getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu);
    virtual void computeStaple(Links *lattice, int i, int j, int k, int l, int mu);
    // Setters
    void setN(int *N);
    void setIndexOrganiser(IndexOrganiser *Index) { m_Index = Index; }
};

#endif // ACTION_H
