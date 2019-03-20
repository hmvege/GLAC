#include "luscheraction.h"
#include "config/parameters.h"
#include "math/functions.h"
#include "parallelization/communicator.h"

LuscherAction::LuscherAction(): Action()
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
}

LuscherAction::~LuscherAction()
{
}

double LuscherAction::getDeltaAction(SU3 U, SU3 UPrime)
{
    return traceRealMultiplication((UPrime - U),m_staple)*m_multiplicationFactor;
}

void LuscherAction::computeStaple(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
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

Lattice<SU3> LuscherAction::getActionDerivative(Lattice<SU3> *lattice, int mu)
{
    /*!
     * A slightly modified version of the Action derivative
     * Multiply the product of this with each of the 8 \f$T^a\f$ generators
     * Take trace, real, and the multiply with \f$\beta/3\f$.
     * Multiply with 8 generators \f$T^a\f$ again, sum and return matrix
     * Can multiply the return matrix exactly instead of doing the 8 matrix multiplications and sums ect.
     * See python script gellmann_matrix_multiplication.py for more details.
     */

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
    
    // Luscher method of taking derivative
    for (unsigned int iSite = 0; iSite < m_latticeStaple.m_latticeSize; iSite++) {
        m_tempStaple2[iSite][0] = 0;
        m_tempStaple2[iSite][1] = (m_latticeStaple[iSite][9] - 2*m_latticeStaple[iSite][1] + m_latticeStaple[iSite][17])*0.3333333333333333;
        m_tempStaple2[iSite][2] =    0.5 * (m_latticeStaple[iSite][6]  - m_latticeStaple[iSite][2]);
        m_tempStaple2[iSite][3] =  - 0.5 * (m_latticeStaple[iSite][3]  + m_latticeStaple[iSite][7]);
        m_tempStaple2[iSite][4] =    0.5 * (m_latticeStaple[iSite][12] - m_latticeStaple[iSite][4]);
        m_tempStaple2[iSite][5] =  - 0.5 * (m_latticeStaple[iSite][5]  + m_latticeStaple[iSite][13]);
        m_tempStaple2[iSite][6] = - m_tempStaple2[iSite][2];
        m_tempStaple2[iSite][7] = m_tempStaple2[iSite][3];
        m_tempStaple2[iSite][8] = 0;
        m_tempStaple2[iSite][9] = (m_latticeStaple[iSite][1] - 2*m_latticeStaple[iSite][9] + m_latticeStaple[iSite][17])*0.3333333333333333;
        m_tempStaple2[iSite][10] =   0.5 * (m_latticeStaple[iSite][14] - m_latticeStaple[iSite][10]);
        m_tempStaple2[iSite][11] = - 0.5 * (m_latticeStaple[iSite][11] + m_latticeStaple[iSite][15]);
        m_tempStaple2[iSite][12] = - m_tempStaple2[iSite][4];
        m_tempStaple2[iSite][13] = m_tempStaple2[iSite][5];
        m_tempStaple2[iSite][14] = - m_tempStaple2[iSite][10];
        m_tempStaple2[iSite][15] = m_tempStaple2[iSite][11];
        m_tempStaple2[iSite][16] = 0;
        m_tempStaple2[iSite][17] = (m_latticeStaple[iSite][1] + m_latticeStaple[iSite][9] - 2*m_latticeStaple[iSite][17])*0.3333333333333333;
    }
    
    return m_tempStaple2;
}
