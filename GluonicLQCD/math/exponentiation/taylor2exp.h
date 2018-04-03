#ifndef TAYLOR2EXP_H
#define TAYLOR2EXP_H

#include "su3exp.h"

class Taylor2Exp : public SU3Exp
{
private:
    SU3 I, QSquared;
public:
    Taylor2Exp();
    SU3 exp(SU3 Q);
};

#endif // TAYLOREXP_H
