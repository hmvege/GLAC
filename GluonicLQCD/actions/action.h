#ifndef ACTION_H
#define ACTION_H

#include "math/Links.h"
#include "parallelization/parallel.h"
#include <vector>

using std::cout;
using std::endl;

class Action
{
protected:
    // Index length array
    unsigned int *m_N;
    // For handling the shift-method in parallelization
    std::vector<int> m_position;
public:
    Action();
    virtual ~Action();
    virtual double getDeltaAction(Links * lattice, SU3 U, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    virtual void computeStaple(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    // Setters
    void setN(unsigned int *N);
    // Getters
    unsigned int * getN() { return m_N; }

    virtual SU3 getActionDerivative(Links * lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
};



#endif // ACTION_H
