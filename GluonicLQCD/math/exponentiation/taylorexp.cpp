#include "taylorexp.h"

TaylorExp::TaylorExp()
{
    I.identity();
}

SU3 TaylorExp::exp(SU3 Q)
{
    /*
     * Exponentiate using regular Taylor expansion.
     */
    QSquared = Q;
    QSquared *= Q;
    QCubed = QSquared;
    QCubed *= Q;
    QQuartic = QCubed;
    QQuartic *= Q;
    return I + Q + QSquared*0.5 + QCubed/6.0 + QQuartic/24.0;
}
