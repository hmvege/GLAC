#include "plaquette.h"
#include "correlator.h"
#include "links.h"
#include "functions.h"
#include "parallelization/indexorganiser.h"
#include <vector>

// TEMP
#include <iostream>
using std::cout;
using std::endl;

Plaquette::Plaquette() : Correlator()
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
    delete [] m_N;
    delete [] muIndex;
    delete [] nuIndex;
}

void Plaquette::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
    multiplicationFactor = 18.0*m_latticeSize;
}

double Plaquette::calculate(Links *lattice)
{
    double gamma = 0;
    SU3 P;
//    for (int i = 0; i < m_N; i++) {
//        for (int j = 0; j < m_N; j++) {
//            for (int k = 0; k < m_N; k++) {
//                for (int l = 0; l < m_N_T; l++) {
//                    for (int mu = 0; mu < 4; mu++) {
//                        for (int nu = mu+1; nu < 4; nu++) {
//                            lorentzIndex(mu,muIndex);
//                            lorentzIndex(nu,nuIndex);
//                            P += lattice[stapleIndex(i,j,k,l,m_N,m_N_T)].U[mu]
//                                    *lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N,m_N_T)].U[nu]
//                                    *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N,m_N_T)].U[mu])
//                                    *inverse(lattice[stapleIndex(i,j,k,l,m_N,m_N_T)].U[nu]);
//                        }
//                    }
//                }
//            }
//        }
//    }
    for (int i = 0; i < m_N[0]; i++) {
        for (int j = 0; j < m_N[1]; j++) {
            for (int k = 0; k < m_N[2]; k++) {
                for (int l = 0; l < m_N[3]; l++) {
                    for (int mu = 0; mu < 4; mu++) {
                        lorentzIndex(mu,muIndex); // Saves quite a few flops by not figuring out the mu index every time
                        for (int nu = mu+1; nu < 4; nu++) {
                            indexes[0] = i;
                            indexes[1] = j;
                            indexes[2] = k;
                            indexes[3] = l;
                            lorentzIndex(nu,nuIndex);
//                            P += lattice[stapleIndex(i,j,k,l,m_N)].U[mu]
//                                    *lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N)].U[nu]
//                                    *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N)].U[mu])
//                                    *inverse(lattice[stapleIndex(i,j,k,l,m_N)].U[nu]);

                            // Need to fix the index back to stapleIndex-like function/class!!

//                            P += lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[mu]
//                                    *lattice[getIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N[1],m_N[2],m_N[3])].U[nu]
//                                    *lattice[getIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N[1],m_N[2],m_N[3])].U[mu].inv()
//                                    *lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[nu].inv();

                            P += lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[mu]
                                    *m_Index->getPositiveLink(lattice,indexes,mu,muIndex,nu)
                                    *m_Index->getPositiveLink(lattice,indexes,nu,nuIndex,mu).inv()
//                                    *lattice[getIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3],m_N[1],m_N[2],m_N[3])].U[nu]
//                                    *lattice[getIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3],m_N[1],m_N[2],m_N[3])].U[mu].inv()
                                    *lattice[getIndex(i,j,k,l,m_N[1],m_N[2],m_N[3])].U[nu].inv();
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
//    return (P.mat[0].re() + P.mat[4].re() + P.mat[8].re())/multiplicationFactor; // 3 from SU3, 6 from number of plaquettes, 3*6=18
    return gamma/multiplicationFactor; // 3 from SU3, 6 from number of plaquettes, 3*6=18
//    return gamma/3.0/6.0/m_latticeSize; // 3 from SU3, 6 from number of plaquettes
}
