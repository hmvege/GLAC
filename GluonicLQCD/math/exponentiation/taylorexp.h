#ifndef TAYLOREXP_H
#define TAYLOREXP_H

#include "su3exp.h"

class TaylorExp : public SU3Exp
{
private:
    // Taylor degree
    unsigned int m_N;
    double m_taylorFactor = 1;
    SU3 m_QMul, m_QSum;
public:
    TaylorExp(unsigned int N);
    SU3 exp(SU3 Q);
};

#endif // TAYLOREXP_H
