#include "plaquette.h"
#include "correlator.h"
#include "links.h"
#include "functions.h"

Plaquette::Plaquette(int N, int N_T) : Correlator(N, N_T)
{
    muIndex = new int[4];
    nuIndex = new int[4];
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
}

Plaquette::~Plaquette()
{
    delete [] muIndex;
    delete [] nuIndex;
}

double Plaquette::calculate(Links *lattice)
{
    double gamma = 0;
    SU3 P;
    // PARALLELIZE HERE
    for (int i = 0; i < m_N; i++) {
        for (int j = 0; j < m_N; j++) {
            for (int k = 0; k < m_N; k++) {
                for (int l = 0; l < m_N_T; l++) {
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = mu+1; nu < 4; nu++) {
                            lorentzIndex(mu,muIndex);
                            lorentzIndex(nu,nuIndex);
                            P += lattice[stapleIndex(i,j,k,l,m_N,m_N_T)].U[mu]
                                    *lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N,m_N_T)].U[nu]
                                    *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N,m_N_T)].U[mu])
                                    *inverse(lattice[stapleIndex(i,j,k,l,m_N,m_N_T)].U[nu]);
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < 3; i++)
    {
        gamma += P.mat[i*3+i].re();
    }
    return gamma/3.0/6.0/m_latticeSize; // 3 from SU3, 6 from number of plaquettes
}
