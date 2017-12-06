#include "wilsongaugeaction.h"
#include "config/parameters.h"
#include "math/functions.h"

//using LatOps::FORWARDS;
//using LatOps::BACKWARDS;
//using LatOps::shift;
//using LatOps::imagTrace;
//using LatOps::subtractImag;
//using LatOps::inv;

WilsonGaugeAction::WilsonGaugeAction(): Action()
{
    m_beta = Parameters::getBeta();
    m_multiplicationFactor = -m_beta/3.0;
    for (int i = 0; i < 4; i++) {
        m_muIndex[i] = 0;
        m_nuIndex[i] = 0;
    }
    m_latticeStaple.allocate(m_N);
    m_tempStaple1.allocate(m_N);
    m_tempStaple2.allocate(m_N);
    m_tempDiag.allocate(m_N);
}

WilsonGaugeAction::~WilsonGaugeAction()
{
}

double WilsonGaugeAction::getDeltaAction(SU3 U, SU3 UPrime)
{
    return traceRealMultiplication((UPrime - U),m_staple)*m_multiplicationFactor;
}

void WilsonGaugeAction::computeStaple(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
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
        m_staple1 = Parallel::Communicator::getPositiveLink(lattice,m_position,mu,m_muIndex,nu);
        m_staple1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,m_nuIndex,mu).inv();
        m_staple1 *= lattice[nu][Parallel::Index::getIndex(i,j,k,l)].inv();
        // Getting second part of staple
        m_staple2 = Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,mu,m_muIndex,nu,m_nuIndex,nu).inv();
        m_staple2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,m_nuIndex,mu).inv();
        m_staple2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,m_nuIndex,nu);
        // Sums staple
        m_staple += m_staple1;
        m_staple += m_staple2;
    }
}

Lattice<SU3> WilsonGaugeAction::getActionDerivative(Lattice<SU3> *lattice, int mu)
{
    // Computes the staple for current link
    m_latticeStaple.zeros();
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        // Retrieves first staple part
        m_tempStaple1 = shift(lattice[nu],FORWARDS,mu);
        // Sets up staple to pass
        m_tempStaple2 = inv(m_tempStaple1);
        m_tempStaple2 *= inv(lattice[mu]);
        m_tempStaple2 *= lattice[nu];
        // Multiplies together local staple
        m_tempStaple1 *= inv(shift(lattice[mu],FORWARDS,nu));
        m_tempStaple1 *= inv(lattice[nu]);
        // Sums the staples
        m_latticeStaple += m_tempStaple1;
        m_latticeStaple += shift(m_tempStaple2,BACKWARDS,nu);
    }

    // Multiplying staple together with link
    m_latticeStaple = lattice[mu]*m_latticeStaple;

    // MORNINGSTAR METHOD
    m_tempStaple1 = inv(m_latticeStaple);
    m_tempStaple1 -= m_latticeStaple;
    m_tempDiag = imagTrace(m_tempStaple1)/3.0;
    m_tempStaple1 = subtractImag(m_tempStaple1,m_tempDiag);
    m_tempStaple1 *= 0.5;
    return m_tempStaple1;
//    Omega = m_staple.inv();
//    Omega -= m_staple;
//    tempDiag = (Omega.mat[1] + Omega.mat[9] + Omega.mat[17])/3.0;
//    for (int i = 1; i < 18; i+=8) { // 8 is subtracting from identity
//        Omega.mat[i] -= tempDiag;
//    }
//    m_X = Omega*0.5;

    // LUSCHER METHOD
    // Multiply the product of this with each of the 8 T^a generators
    // Take trace, real, and the multiply with beta/3
    // Multiply with 8 generators T^a again, sum and return matrix
    // Can multiply the return matrix exactly instead of doing the 8 matrix multiplications and sums ect. See python script gellmann_matrix_multiplication.py.
//    // Diagonals
//    m_X.mat[0] = 0;
//    m_X.mat[1] = (m_staple.mat[9] - 2*m_staple.mat[1] + m_staple.mat[17])*0.3333333333333333;
//    m_X.mat[8] = 0;
//    m_X.mat[9] = (m_staple.mat[1] - 2*m_staple.mat[9] + m_staple.mat[17])*0.3333333333333333;
//    m_X.mat[16] = 0;
//    m_X.mat[17] = (m_staple.mat[1] + m_staple.mat[9] - 2*m_staple.mat[17])*0.3333333333333333;
//    // Upper triangular
//    m_X.mat[2] =    0.5 * (m_staple.mat[6]  - m_staple.mat[2]);
//    m_X.mat[3] =  - 0.5 * (m_staple.mat[3]  + m_staple.mat[7]);
//    m_X.mat[4] =    0.5 * (m_staple.mat[12] - m_staple.mat[4]);
//    m_X.mat[5] =  - 0.5 * (m_staple.mat[5]  + m_staple.mat[13]);
//    m_X.mat[10] =   0.5 * (m_staple.mat[14] - m_staple.mat[10]);
//    m_X.mat[11] = - 0.5 * (m_staple.mat[11] + m_staple.mat[15]);
//    // Lower triangular
//    m_X.mat[6] = - m_X.mat[2];
//    m_X.mat[7] = m_X.mat[3];
//    m_X.mat[12] = - m_X.mat[4];
//    m_X.mat[13] = m_X.mat[5];
//    m_X.mat[14] = - m_X.mat[10];
//    m_X.mat[15] = m_X.mat[11];

//    return m_X;
}
