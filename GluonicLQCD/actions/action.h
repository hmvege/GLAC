#ifndef ACTION_H
#define ACTION_H

#include "math/lattice.h"
#include <vector>

using std::cout;
using std::endl;

class Action
{
protected:
    // Index length array
    std::vector<unsigned int> m_N;
    // For handling the shift-method in parallelization
    std::vector<int> m_position;
public:
    Action();
    virtual ~Action();
    virtual double getDeltaAction(SU3 U, SU3 UPrime);
    virtual void computeStaple(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    virtual Lattice<SU3> getActionDerivative(Lattice<SU3> * lattice, int mu);
};



#endif // ACTION_H
