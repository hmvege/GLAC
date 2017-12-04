#ifndef ACTION_H
#define ACTION_H

#include "math/lattice.h"
#include "parallelization/parallel.h"
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
    // Setters
    void setN(std::vector<unsigned int> N);
    // Getters
    std::vector<unsigned int> getN() { return m_N; }

    virtual Lattice<SU3> getActionDerivative(Lattice<SU3> * lattice, int mu);
};



#endif // ACTION_H
