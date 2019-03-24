/*!
 * \class TaylorExp
 *
 * \brief uses a general Taylor polynomial to approximate the exponentiation of a SU(3) matrix.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef TAYLOREXP_H
#define TAYLOREXP_H

#include "su3exp.h"

class TaylorExp : public SU3Exp
{
private:
    // Taylor degree
    unsigned int m_N;
    double m_taylorFactor = 1;
    SU3 m_QMul, m_QSum;
public:
    TaylorExp(unsigned int N);
    SU3 exp(SU3 Q);

    void setTaylorDegree(unsigned int N);
};

#endif // TAYLOREXP_H
