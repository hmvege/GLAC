#include "wilsongaugeaction.h"
#include "action.h"
#include "functions.h"

WilsonGaugeAction::WilsonGaugeAction(int N, double beta): Action(N)
{
    m_beta = beta;
    muIndex = new int[4];
    nuIndex = new int[4];
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
}

WilsonGaugeAction::~WilsonGaugeAction()
{
    delete [] muIndex;
    delete [] nuIndex;
}

double WilsonGaugeAction::getDeltaAction(Links *lattice, SU3 U, int i, int j, int k, int l, int mu)
{
    double S = 0;
    SU3 A;
    SU3 tr;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        lorentzIndex(mu,muIndex);
        lorentzIndex(nu,nuIndex);

        A += lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], m_N)].U[nu]
                *inverse(lattice[stapleIndex(i+muIndex[0]+nuIndex[0],j+muIndex[1]+nuIndex[1],k+muIndex[2]+nuIndex[2],l+muIndex[3]+nuIndex[3], m_N)].U[mu])
                *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3], m_N)].U[nu])
                + inverse(lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], m_N)].U[nu])
                *inverse(lattice[stapleIndex(i+muIndex[0]-nuIndex[0],j+muIndex[1]-nuIndex[1],k+muIndex[2]-nuIndex[2],l+muIndex[3]-nuIndex[3], m_N)].U[mu])
                *lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3], m_N)].U[nu];
    }
    tr = (U - lattice[index(i,j,k,l,m_N)].U[mu])*A;
    for (int i = 0; i < 3; i++)
    {
        S += tr.mat[i*3+i].re;
    }
    return m_beta/3.0*S; // Should be N=3 as in the Gauge symmetry?
}
