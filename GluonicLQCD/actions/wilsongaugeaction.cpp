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
    m_indexes[0] = i;
    m_indexes[1] = j;
    m_indexes[2] = k;
    m_indexes[3] = l;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        updateNuIndex(nu);

        m_staple += m_Index->getPositiveLink(lattice,m_indexes,mu,m_muIndex,nu)
                *m_Index->getPositiveLink(lattice,m_indexes,nu,m_nuIndex,mu).inv()
                *lattice[m_Index->getIndex(i,j,k,l)].U[nu].inv()
                + m_Index->getNeighboursNeighbourLink(lattice,m_indexes,mu,m_muIndex,nu,m_nuIndex,nu).inv()
                *m_Index->getNegativeLink(lattice,m_indexes,nu,m_nuIndex,mu).inv()
                *m_Index->getNegativeLink(lattice,m_indexes,nu,m_nuIndex,nu);
    }
}

SU3 WilsonGaugeAction::getActionDerivative(Links * lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    m_staple.zeros();
//    m_X.zeros();
    updateMuIndex(mu);
    m_indexes[0] = i;
    m_indexes[1] = j;
    m_indexes[2] = k;
    m_indexes[3] = l;
    // First: calculate the X_mu(see notes)
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        updateNuIndex(nu);
        m_staple += m_Index->getPositiveLink(lattice,m_indexes,mu,m_muIndex,nu)
                *m_Index->getPositiveLink(lattice,m_indexes,nu,m_nuIndex,mu).inv()
                *lattice[m_Index->getIndex(i,j,k,l)].U[nu].inv()
                + m_Index->getNeighboursNeighbourLink(lattice,m_indexes,mu,m_muIndex,nu,m_nuIndex,nu).inv()
                *m_Index->getNegativeLink(lattice,m_indexes,nu,m_nuIndex,mu).inv()
                *m_Index->getNegativeLink(lattice,m_indexes,nu,m_nuIndex,nu);
    }
    // Multiply the product of this with each of the 8 T^a generators
    // Take trace, real, and the multiply with beta/3
    // Multiply with 8 generators T^a again, sum and return matrix

//    SU3 m_staplez = lattice[m_Index->getIndex(i,j,k,l)].U[mu]*m_staple;
//    SU3 Q, C;
//    double temp;
//    temp = (m_staplez.inv() - m_staplez).trace().im()/6.0;
//    C.zeros();
//    C.mat[1] = temp;
//    C.mat[9] = temp;
//    C.mat[17] = temp;

//    Q = (m_staplez.inv() - m_staplez)*0.5 - C;

//    m_X = Q*0.5;
//    m_X.print();

    // Multiply X_mu with the U_mu
    m_staple = lattice[m_Index->getIndex(i,j,k,l)].U[mu] * m_staple;


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

//    m_X.printMachine();
//    m_X.makeHermitian(); // Matrix now Hermitian

//    std::cout << "printing X" << std::endl;
//    m_X.printMachine();
//    std::cout << std::endl;
//    m_X.makeHermitian(); // Matrix now anti-Hermitian
//    m_X.printMachine();
//    exit(1);
    // No multiplication factor needed as the -3/beta is cancelled in the Z(W_i).
    return m_X;//*m_multiplicationFactor;
}
