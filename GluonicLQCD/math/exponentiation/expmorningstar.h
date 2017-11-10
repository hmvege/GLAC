#ifndef EXPMORNINGSTAR_H
#define EXPMORNINGSTAR_H

#include "su3exp.h"

class ExpMorningstar : public SU3Exp
{
private:
    // Base constants
    complex f[3];
    complex h[3];
    double c0, c1, u, w, theta, xi0, c0max;
    SU3 f0; // Matrix for filling return matrix
    // Derived constants
    SU3 QSquared, QCubed;
    double cosu, cosw, sinu, sinw, uu, ww, sin2u, cos2u;
public:
    ExpMorningstar();
    SU3 exp(SU3 Q);
};

#endif // EXPMORNINGSTAR_H
