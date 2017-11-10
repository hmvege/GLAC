#ifndef TAYLOREXP_H
#define TAYLOREXP_H

#include "su3exp.h"

class TaylorExp : public SU3Exp
{
private:
    SU3 I, QSquared, QCubed, QQuartic;
public:
    TaylorExp();
    SU3 exp(SU3 Q);
};

#endif // TAYLOREXP_H
