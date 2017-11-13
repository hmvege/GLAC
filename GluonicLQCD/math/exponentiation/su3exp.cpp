#include "su3exp.h"

SU3Exp::SU3Exp()
{
    /*
     * Default SU3 exponentiation is with the Morningstar method.
     */
}

SU3Exp::~SU3Exp()
{

}

SU3 SU3Exp::exp(SU3 Q)
{
    /*
     * Takes the exponential of an Hermitian matrix using the Morningstar method.
     */
    f0.zeros();
    Q.makeHermitian();

    QSquared = Q*Q;
    QCubed = Q*QSquared;

    c0 = 0.3333333333333333*(QCubed.mat[0] + QCubed.mat[8] + QCubed.mat[16]);
    c1 = 0.5*(QSquared.mat[0] + QSquared.mat[8] + QSquared.mat[16]);

    c0max = 0.6666666666666666 * c1 * sqrt(c1 * 0.3333333333333333);
    theta = acos(fabs(c0)/c0max);
    u = sqrt(0.3333333333333333*c1) * cos(0.3333333333333333 * theta);
    w = sqrt(c1) * sin(0.3333333333333333 * theta);

    // Sets shortenings
    uu = u*u;
    ww = w*w;
    cosu = cos(u);
    sinu = sin(u);
    cosw = cos(w);
    sinw = sin(w);
    cos2u = cos(2*u);
    sin2u = sin(2*u);

    // Find xi(w)
    if (fabs(w) > 0.05) {
        xi0 = sinw/w;
    } else {
        xi0 = 1.0 - ww*(1 - 0.05*ww*(1 - ww/42.0))/6.0; // 1/6.0 and 1/42.0
    }

    // Sets h
    h[0].setRe(8*uu*cosu*cosw + 2*u*xi0*(3*uu + ww)*sinu + (uu - ww)*cos2u);
    h[0].setIm(-8*uu*sinu*cosw + 2*u*xi0*(3*uu + ww)*cosu + (uu - ww)*sin2u);
    h[1].setRe(-2*u*cosu*cosw + 2*u*cos2u + xi0*(3*uu - ww)*sinu);
    h[1].setIm(2*u*sinu*cosw + 2*u*sin2u + xi0*(3*uu - ww)*cosu);
    h[2].setRe(-3*u*xi0*sinu - cosu*cosw + cos2u);
    h[2].setIm(-3*u*xi0*cosu + sinu*cosw + sin2u);

    // Sets f
    for (int i = 0; i < 3; i++) {
        f[i] = h[i] / (9*uu - ww);
    }

    // Checks for negative c0 coefficient
    if (c0 < 0) {
        f[0].conjugate();
        f[1] = -f[1].c();
        f[2].conjugate();
    }

    // Sets the first matrix, I*f0
    f0.setComplex(f[0],0);
    f0.setComplex(f[0],8);
    f0.setComplex(f[0],16);

    return f0 + Q*f[1] + QSquared*f[2];
}
