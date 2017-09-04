#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "links.h"

class Correlator
{
protected:
//    int m_N;
//    int m_N_T;
    int *m_N;
    double m_latticeSize;
public:
//    Correlator(int N, int N_T);
    Correlator();
    ~Correlator();
    virtual double calculate(Links *lattice);
    void setN(int *N);
    void setLatticeSize(int latticeSize) { m_latticeSize = double(latticeSize); }
};

#endif // CORRELATOR_H
