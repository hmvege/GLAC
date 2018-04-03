#include "taylor2exp.h"

Taylor2Exp::Taylor2Exp()
{
    I.identity();
}

SU3 Taylor2Exp::exp(SU3 Q)
{
    /*
     * Exponentiate using regular Taylor expansion to the second order.
     */
    QSquared = Q;
    QSquared *= Q;
    return I + Q + QSquared*0.5;
}
