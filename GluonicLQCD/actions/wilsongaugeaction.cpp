#include "wilsongaugeaction.h"
#include "action.h"
#include "functions.h"
#include <vector>

WilsonGaugeAction::WilsonGaugeAction(double beta): Action()
{
    m_beta = beta;
    multiplicationFactor = -m_beta/3.0;
//    muIndex = new int[4];
//    nuIndex = new int[4];
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
}

WilsonGaugeAction::~WilsonGaugeAction()
{
//    delete [] muIndex;
//    delete [] nuIndex;
}

double WilsonGaugeAction::getDeltaAction(Links *lattice, SU3 UPrime, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
//    tr = (UPrime - lattice[m_Index->getIndex(i,j,k,l)].U[mu])*staple;
//    return (tr.mat[0] + tr.mat[8] + tr.mat[16])*multiplicationFactor; // Should be N=3 as in the Gauge symmetry
    return traceRealMultiplication((UPrime - lattice[m_Index->getIndex(i,j,k,l)].U[mu]),staple)*multiplicationFactor;
}

void WilsonGaugeAction::computeStaple(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    staple.zeros();
//    lorentzIndex(mu,muIndex);
    updateMuIndex(mu);
    indexes[0] = i;
    indexes[1] = j;
    indexes[2] = k;
    indexes[3] = l;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
//        lorentzIndex(nu,nuIndex);
        updateNuIndex(nu);

        // Need to fix the index back to stapleIndex-like function/class!!
        staple += m_Index->getPositiveLink(lattice,indexes,mu,muIndex,nu)
                *m_Index->getPositiveLink(lattice,indexes,nu,nuIndex,mu).inv()
                *lattice[m_Index->getIndex(i,j,k,l)].U[nu].inv()
                + m_Index->getNeighboursNeighbourLink(lattice,indexes,mu,muIndex,nu,nuIndex,nu).inv()
                *m_Index->getNegativeLink(lattice,indexes,nu,nuIndex,mu).inv()
                *m_Index->getNegativeLink(lattice,indexes,nu,nuIndex,nu);
    }
}
