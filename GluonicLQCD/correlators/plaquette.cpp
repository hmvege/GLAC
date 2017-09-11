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
    setLatticeSize(N*N*N*N_T);
}

Plaquette::~Plaquette()
{
    delete [] muIndex;
    delete [] nuIndex;
}

void Plaquette::setLatticeSize(int latticeSize) {
    m_latticeSize = latticeSize;
    m_multiplicationFactor = 1/(3.0*6.0*double(m_latticeSize));
}

double Plaquette::calculate(Links *lattice)
{
    double gamma = 0;
    SU3 P;
    for (int i = 0; i < m_N; i++) {
        for (int j = 0; j < m_N; j++) {
            for (int k = 0; k < m_N; k++) {
                for (int l = 0; l < m_N_T; l++) {
                    for (int mu = 0; mu < 4; mu++) {
                        lorentzIndex(mu,muIndex);
                        for (int nu = mu+1; nu < 4; nu++) {
                            lorentzIndex(nu,nuIndex);
                            P += lattice[stapleIndex(i,j,k,l,m_N,m_N_T)].U[mu]
                                    *lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N,m_N_T)].U[nu]
                                    *lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N,m_N_T)].U[mu].inv()
                                    *lattice[stapleIndex(i,j,k,l,m_N,m_N_T)].U[nu].inv();
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < 3; i++)
    {
        gamma += P.mat[i*3+i].re;
    }
    return gamma*m_multiplicationFactor; // 3 from SU3, 6 from number of plaquettes
}
