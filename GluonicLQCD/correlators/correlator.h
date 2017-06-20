#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "links.h"

class Correlator
{
protected:
    int m_N;
    int m_latticeSize;
public:
    Correlator(int N);
    ~Correlator();
    virtual double calculate(Links *lattice);
};

#endif // CORRELATOR_H
