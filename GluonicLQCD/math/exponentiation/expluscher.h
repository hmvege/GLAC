#ifndef EXPLUSCHER_H
#define EXPLUSCHER_H

#include "su3exp.h"

class ExpLuscher : public SU3Exp
{
private:
    SU3 U1, U2, U3, E;
    complex x1, x2, x3, div1, div2, div3, X1221X, X1331X, X2332X, sqrdFactor;
public:
    ExpLuscher();
    SU3 exp(SU3 Q);
};

#endif // EXPLUSCHER_H
