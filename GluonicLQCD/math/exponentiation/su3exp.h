#ifndef SU3EXP_H
#define SU3EXP_H

#include "math/matrices/su3.h"

class SU3Exp
{
public:
    SU3Exp();
    virtual ~SU3Exp();
    virtual SU3 exp(SU3 Q);
};

#endif // SU3EXP_H
