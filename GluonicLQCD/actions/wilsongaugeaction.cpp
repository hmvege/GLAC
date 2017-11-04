#include "wilsongaugeaction.h"
#include "action.h"
#include "functions.h"
#include <vector>

WilsonGaugeAction::WilsonGaugeAction(double beta): Action()
{
    m_beta = beta;
    m_multiplicationFactor = -m_beta/3.0;
    for (int i = 0; i < 4; i++) {
        m_muIndex[i] = 0;
        m_nuIndex[i] = 0;
    }
}

WilsonGaugeAction::~WilsonGaugeAction()
{
}

double WilsonGaugeAction::getDeltaAction(Links *lattice, SU3 UPrime, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    return traceRealMultiplication((UPrime - lattice[m_Index->getIndex(i,j,k,l)].U[mu]),m_staple)*m_multiplicationFactor;
}

void WilsonGaugeAction::computeStaple(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    m_staple.zeros();
    updateMuIndex(mu);
    m_position[0] = i;
    m_position[1] = j;
    m_position[2] = k;
    m_position[3] = l;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        updateNuIndex(nu);
        // Getting first part of staple
        m_staple1 = m_Index->getPositiveLink(lattice,m_position,mu,m_muIndex,nu);
        m_staple1 *= m_Index->getPositiveLink(lattice,m_position,nu,m_nuIndex,mu).inv();
        m_staple1 *= lattice[m_Index->getIndex(i,j,k,l)].U[nu].inv();
        // Getting second part of staple
        m_staple2 = m_Index->getNeighboursNeighbourLink(lattice,m_position,mu,m_muIndex,nu,m_nuIndex,nu).inv();
        m_staple2 *= m_Index->getNegativeLink(lattice,m_position,nu,m_nuIndex,mu).inv();
        m_staple2 *= m_Index->getNegativeLink(lattice,m_position,nu,m_nuIndex,nu);
        // Sums staple
        m_staple += m_staple1;
        m_staple += m_staple2;
    }
}

SU3 WilsonGaugeAction::getActionDerivative(Links * lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    computeStaple(lattice,i,j,k,l,mu);
    m_staple = lattice[m_Index->getIndex(i,j,k,l)].U[mu]*m_staple;

    // MORNINGSTAR METHOD
//    Omega = m_staple.inv();
//    Omega -= m_staple;
//    tempDiag = Omega.trace().im()/6.0;
//    tempDiag = (Omega.mat[1] + Omega.mat[9] + Omega.mat[17])/6.0;
//    Q = Omega*0.5;
//    Q.mat[1] -= tempDiag;
//    Q.mat[9] -= tempDiag;
//    Q.mat[17] -= tempDiag;
//    m_X = Q*0.5;

    // LUSCHER METHOD
    // Multiply the product of this with each of the 8 T^a generators
    // Take trace, real, and the multiply with beta/3
    // Multiply with 8 generators T^a again, sum and return matrix
    // Can multiply the return matrix exactly instead of doing the 8 matrix multiplications and sums ect. See python script gellmann_matrix_multiplication.py.
    // Diagonals
    m_X.mat[0] = 0;
    m_X.mat[1] = (m_staple.mat[9] - 2*m_staple.mat[1] + m_staple.mat[17])*0.1666666666666666;
    m_X.mat[8] = 0;
    m_X.mat[9] = (m_staple.mat[1] - 2*m_staple.mat[9] + m_staple.mat[17])*0.1666666666666666;
    m_X.mat[16] = 0;
    m_X.mat[17] = (m_staple.mat[1] + m_staple.mat[9] - 2*m_staple.mat[17])*0.1666666666666666;

    // Upper triangular
    m_X.mat[2] =    0.25 * (m_staple.mat[6]  - m_staple.mat[2]);
    m_X.mat[3] =  - 0.25 * (m_staple.mat[3]  + m_staple.mat[7]);
    m_X.mat[4] =    0.25 * (m_staple.mat[12] - m_staple.mat[4]);
    m_X.mat[5] =  - 0.25 * (m_staple.mat[5]  + m_staple.mat[13]);
    m_X.mat[10] =   0.25 * (m_staple.mat[14] - m_staple.mat[10]);
    m_X.mat[11] = - 0.25 * (m_staple.mat[11] + m_staple.mat[15]);

    // Lower triangular
    m_X.mat[6] = - m_X.mat[2];
    m_X.mat[7] = m_X.mat[3];
    m_X.mat[12] = - m_X.mat[4];
    m_X.mat[13] = m_X.mat[5];
    m_X.mat[14] = - m_X.mat[10];
    m_X.mat[15] = m_X.mat[11];

    return m_X;;
}
