/*!
 * \class SU3Exp
 *
 * \brief SU3 matrix exponentiation method that uses the method from https://journals.aps.org/prd/abstract/10.1103/PhysRevD.69.054501
 *
 * Also serves as base class for the other exponentiation methods.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef SU3EXP_H
#define SU3EXP_H

#include "math/matrices/su3.h"

class SU3Exp
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
    SU3Exp();
    virtual ~SU3Exp();
    virtual SU3 exp(SU3 Q);
};

#endif // SU3EXP_H
