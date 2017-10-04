#include "plaquette.h"
#include <vector>
#include "correlator.h"
#include "links.h"
#include "functions.h"
#include "parallelization/indexorganiser.h"

Plaquette::Plaquette() : Correlator()
{
//    muIndex = new int[4];
//    nuIndex = new int[4];
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
}

Plaquette::~Plaquette()
{
//    delete [] muIndex;
//    delete [] nuIndex;
}

void Plaquette::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
    multiplicationFactor = 18.0*m_latticeSize;
}

double Plaquette::calculate(Links *lattice)
{
    P.zeros();
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    indexes[0] = i;
                    indexes[1] = j;
                    indexes[2] = k;
                    indexes[3] = l;
                    for (int mu = 0; mu < 4; mu++) {
//                        lorentzIndex(mu,muIndex); // Saves quite a few flops by not figuring out the mu index every time
                        updateMuIndex(mu); // Inline function
                        for (int nu = mu+1; nu < 4; nu++) {
//                            lorentzIndex(nu,nuIndex);
                            updateNuIndex(nu); // Inline function
                            P += lattice[m_Index->getIndex(i,j,k,l)].U[mu]
                                    *m_Index->getPositiveLink(lattice,indexes,mu,muIndex,nu)
                                    *m_Index->getPositiveLink(lattice,indexes,nu,nuIndex,mu).inv()
                                    *lattice[m_Index->getIndex(i,j,k,l)].U[nu].inv();
                        }
                    }
                }
            }
        }
    }
    // Small redo here!
    return (P.mat[0].re() + P.mat[4].re() + P.mat[8].re())/multiplicationFactor; // 3 from SU3, 6 from number of plaquettes, 3*6=18
}
