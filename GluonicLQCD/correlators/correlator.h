#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "links.h"

class Correlator
{
protected:
    int m_N;
    int m_N_T;
    int m_latticeSize;
public:
    Correlator(int N, int N_T);
    ~Correlator();
    virtual double calculate(Links *lattice);
    virtual void setLatticeSize(int latticeSize) { m_latticeSize = latticeSize; }
};

#endif // CORRELATOR_H
