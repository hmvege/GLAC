#include "expluscher.h"

ExpLuscher::ExpLuscher()
{

}

/*!
 * \brief ExpLuscher::exp method for exponentiation, found in paper appendix A, https://arxiv.org/abs/hep-lat/0409106
 * \param Q matrix to exponentiate.
 * \return The exponentiated SU3 Q matrix.
 */
SU3 ExpLuscher::exp(SU3 Q)
{
    /*
     * Exponentiation using the Luscher method.
     */
    // Ensures the U's are at zero for later filling
    U1.zeros();
    U2.zeros();
    U3.zeros();

    // Sets elements of the Y matrices
    x1 = (Q.get(0,0) - Q.get(1,1))*0.3333333333333333;
    x2 = (Q.get(0,0) - Q.get(2,2))*0.3333333333333333;
    x3 = (Q.get(1,1) - Q.get(2,2))*0.3333333333333333;

    // Sets often used factors
    X1221X = Q.get(0,1)*Q.get(1,0);
    X1331X = Q.get(0,2)*Q.get(2,0);
    X2332X = Q.get(1,2)*Q.get(2,1);

    // Sets complex denominators
    div1 = (X1221X + x1*x1)*0.0625 - 1.0;
    div2 = (X1331X + x2*x2)*0.0625 - 1.0;
    div3 = (X2332X + x3*x3)*0.25 - 1.0;

    // Setting U1 matrix
    sqrdFactor = x1*0.25 - 1.0;
    U1.setComplex((X1221X*0.0625 + x1*x1*0.0625 + x1*0.5 + 1.0)/div1*(-1.0),0);
    U1.setComplex((Q.get(0,1)*(-0.5))/div1,2);
    U1.setComplex((Q.get(1,0)*(-0.5))/div1,6);
    U1.setComplex((X1221X*0.0625 + sqrdFactor*sqrdFactor)/div1*(-1.0),8);
    U1.mat[16] = 1.0;

    // Setting U2 matrix
    sqrdFactor = x2*0.25 - 1.0;
    U2.setComplex((X1331X*0.0625 + x2*x2*0.0625 + x2*0.5 + 1.0)/div2*(-1),0);
    U2.setComplex((Q.get(0,2)*(-0.5))/div2,4);
    U2.mat[8] = 1.0;
    U2.setComplex((Q.get(2,0)*(-0.5))/div2,12);
    U2.setComplex((X1331X*0.0625 + sqrdFactor*sqrdFactor)/div2*(-1.0),16);

    // Setting U3 matrix
    sqrdFactor = x3*0.5 - 1.0;
    U3.mat[0] = 1.0;
    U3.setComplex((X2332X*0.25 + x3*x3*0.25 + x3 + 1.0)/div3*(-1.0),8);
    U3.setComplex(Q.get(1,2)/div3*(-1.0),10);
    U3.setComplex(Q.get(2,1)/div3*(-1.0),14);
    U3.setComplex((X2332X*0.25 + sqrdFactor*sqrdFactor)/div3*(-1.0),16);

    E = U1;
    E *= U2;
    E *= U3;
    E *= U2;
    E *= U1;
    return E;
}
