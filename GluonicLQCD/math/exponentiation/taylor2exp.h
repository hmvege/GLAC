/*!
 * \class Taylor2Exp
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
