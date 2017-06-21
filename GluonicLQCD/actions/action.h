#ifndef ACTION_H
#define ACTION_H

#include "links.h"
#include "functions.h"

class Action
{
protected:
//    double a;
    int m_N;
    // Lorentz indices arrays
    int *muIndex;
    int *nuIndex;
public:
    Action();
    Action(int latticeSize);
    virtual ~Action();
    virtual double getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu);
    virtual void computeStaple(Links *lattice, int i, int j, int k, int l, int mu);
    // Setters
//    void setLatticeSpacing(double new_a) { a = new_a; }
    void setNLatticePoints(int N) { m_N = N; }
    // Getters
//    double getLatticeSpacing() { return a; }
    int getNLatticePoints() { return m_N; }
};

#endif // ACTION_H
