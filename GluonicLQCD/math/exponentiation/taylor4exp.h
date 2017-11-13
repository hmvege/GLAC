#ifndef TAYLOR4EXP_H
#define TAYLOR4EXP_H

#include "su3exp.h"

class Taylor4Exp : public SU3Exp
{
private:
    SU3 I, QSquared, QCubed, QQuartic;
public:
    Taylor4Exp();
    SU3 exp(SU3 Q);
};

#endif // TAYLOR4EXP_H
