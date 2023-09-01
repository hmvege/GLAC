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

    /** 
     * Method which verifies the Taylor degree
     * 
     * Minimum degree is 1, and maximum is 16 (however, after degree 12,
     * no further improvement are seen)
     */
    void verifyTaylorDegree(int n);
    
public:
    TaylorExp(unsigned int N);

    /** \brief Exponentiate using regular Taylor expansion */
    SU3 exp(SU3 Q);

    void setTaylorDegree(unsigned int N);
};

#endif // TAYLOREXP_H
