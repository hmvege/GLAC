#include "wilsongaugeaction.h"
#include "action.h"
#include "functions.h"

#include <vector>

// TEMP
#include <cmath>
#include <iostream>

WilsonGaugeAction::WilsonGaugeAction(double beta): Action()
{
    m_beta = beta;
    multiplicationFactor = -m_beta/3.0;
    muIndex = new int[4];
    nuIndex = new int[4];
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
}

WilsonGaugeAction::~WilsonGaugeAction()
{
    delete [] m_N;
    delete [] muIndex;
    delete [] nuIndex;
}

double WilsonGaugeAction::getDeltaAction(Links *lattice, SU3 UPrime, int i, int j, int k, int l, int mu)
{
    double S = 0;
    SU3 tr;
    tr = (UPrime - lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[mu])*m_staple;
    for (int n = 0; n < 3; n++)
    {
        S += tr.mat[n*3+n].re();
    }
    return S*multiplicationFactor; // Should be N=3 as in the Gauge symmetry
}

void WilsonGaugeAction::computeStaple(Links *lattice, int i, int j, int k, int l, int mu)
{
    m_staple.zeros();
    lorentzIndex(mu,muIndex);
    indexes[0] = i;
    indexes[1] = j;
    indexes[2] = k;
    indexes[3] = l;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        lorentzIndex(nu,nuIndex);

        // Need to fix the index back to stapleIndex-like function/class!!
        m_staple += m_Index->getPositiveLink(lattice,indexes,mu,muIndex,nu)
                *m_Index->getPositiveLink(lattice,indexes,nu,nuIndex,mu).inv()
                *lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[nu].inv()
                + m_Index->getNeighboursNeighbourLink(lattice,indexes,mu,muIndex,nu,nuIndex,nu).inv()
                *m_Index->getNegativeLink(lattice,indexes,mu,nuIndex,mu).inv()
                *m_Index->getNegativeLink(lattice,indexes,mu,muIndex,mu);

//        m_staple += lattice[getIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N[1],m_N[2],m_N[3])].U[nu] // Gattinger using verbatim staple definition OLD METHOD
//                *lattice[getIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N[1],m_N[2],m_N[3])].U[mu].inv()
//                *lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[nu].inv()
//                + lattice[getIndex(i+muIndex[0]-nuIndex[0],j+muIndex[1]-nuIndex[1],k+muIndex[2]-nuIndex[2],l+muIndex[3]-nuIndex[3],m_N[1],m_N[2],m_N[3])].U[nu].inv()
//                *lattice[getIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3],m_N[1],m_N[2],m_N[3])].U[mu].inv()
//                *lattice[getIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3],m_N[1],m_N[2],m_N[3])].U[nu];

//        m_staple += lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N)].U[nu] // Gattinger using verbatim staple definition OLD METHOD
//                *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N)].U[mu])
//                *inverse(lattice[stapleIndex(i,j,k,l,m_N)].U[nu])
//                + inverse(lattice[stapleIndex(i+muIndex[0]-nuIndex[0],j+muIndex[1]-nuIndex[1],k+muIndex[2]-nuIndex[2],l+muIndex[3]-nuIndex[3],m_N)].U[nu])
//                *inverse(lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3],m_N)].U[mu])
//                *lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3],m_N)].U[nu];
    }
}
