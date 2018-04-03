#include "wilsongaugeaction.h"
#include "config/parameters.h"
#include "math/functions.h"

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
}
