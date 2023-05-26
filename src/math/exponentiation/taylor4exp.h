/*!
 * \class Taylor4Exp
 *
 * \brief Taylor expansion of second order.
 *
 * \todo remove this method, and instead use taylorExp.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
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
