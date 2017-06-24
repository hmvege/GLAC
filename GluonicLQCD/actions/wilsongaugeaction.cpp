#include "wilsongaugeaction.h"
#include "action.h"
#include "functions.h"

WilsonGaugeAction::WilsonGaugeAction(int N, int N_T, double beta): Action(N, N_T)
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

double WilsonGaugeAction::getDeltaAction(Links *lattice, SU3 UPrime, int i, int j, int k, int l, int mu)
{
    double S = 0;
    SU3 tr;
//    m_staple.print();
//    std::cout << i << " "<<j <<" "<< k <<" "<< l <<" " <<mu<< std::endl;
    tr = (UPrime - lattice[index(i,j,k,l,m_N)].U[mu])*m_staple;
    for (int i = 0; i < 3; i++)
    {
        S += tr.mat[i*3+i].re;
    }
    return -m_beta/3.0*S; // Should be N=3 as in the Gauge symmetry?
}

void WilsonGaugeAction::computeStaple(Links *lattice, int i, int j, int k, int l, int mu)
{
    m_staple.zeros();
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
        lorentzIndex(mu,muIndex);
        lorentzIndex(nu,nuIndex);
        m_staple += lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], m_N, m_N_T)].U[nu] // Gattinger using verbatim staple definition
                *inverse(lattice[stapleIndex(i+muIndex[0]+nuIndex[0],j+muIndex[1]+nuIndex[1],k+muIndex[2]+nuIndex[2],l+muIndex[3]+nuIndex[3], m_N, m_N_T)].U[mu])
                *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3], m_N, m_N_T)].U[nu])
                + inverse(lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], m_N, m_N_T)].U[nu])
                *inverse(lattice[stapleIndex(i+muIndex[0]-nuIndex[0],j+muIndex[1]-nuIndex[1],k+muIndex[2]-nuIndex[2],l+muIndex[3]-nuIndex[3], m_N, m_N_T)].U[mu])
                *lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3], m_N, m_N_T)].U[nu];
//        m_staple += lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], m_N, m_N_T)].U[nu] // Gattinger using the link definition
//                *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3], m_N, m_N_T)].U[mu])
//                *inverse(lattice[stapleIndex(i,j,k,l, m_N, m_N_T)].U[nu])
//                + inverse(lattice[stapleIndex(i+muIndex[0]-nuIndex[0],j+muIndex[1]-nuIndex[1],k+muIndex[2]-nuIndex[2],l+muIndex[3]-nuIndex[3], m_N, m_N_T)].U[nu])
//                *inverse(lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3], m_N, m_N_T)].U[mu])
//                *lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3], m_N, m_N_T)].U[nu];
    }
//    for (int i = 0; i < 9; i++) { m_staple.mat[i] /= 6; }
//    staple.print();
//    exit(1);
}
