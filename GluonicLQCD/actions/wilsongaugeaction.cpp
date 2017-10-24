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
    // Multiply X_mu with the U_mu
    m_X = lattice[m_Index->getIndex(i,j,k,l)].U[mu] * m_staple;

    // Multiply the product of this with each of the 8 T^a generators

    // Take trace, real, and the multiply with beta/3

    // Multiply with 8 generators T^a again, sum and return matrix

    return m_X*m_multiplicationFactor;
}
