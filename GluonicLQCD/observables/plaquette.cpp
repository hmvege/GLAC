#include "plaquette.h"
#include <vector>
#include "math/links.h"
#include "math/functions.h"
#include "parallelization/index.h"

Plaquette::Plaquette() : Correlator()
{
}

Plaquette::~Plaquette()
{
}

void Plaquette::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
    m_multiplicationFactor = 18.0*m_latticeSize; // 3 from SU3, 6 from number of plaquettes, 3*6=18
}

double Plaquette::calculate(Links *lattice)
{
    P.zeros();
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    for (int mu = 0; mu < 4; mu++) {
                        updateMuIndex(mu); // Inline function
                        for (int nu = mu+1; nu < 4; nu++) {
                            updateNuIndex(nu); // Inline function
                            PTemp = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
                            PTemp *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
                            PTemp *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
                            PTemp *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
                            P += PTemp;
                        }
                    }
                }
            }
        }
    }
    return (P.mat[0] + P.mat[8] + P.mat[16])/m_multiplicationFactor;
}

double Plaquette::calculate(SU3 *plaquetteStaples)
{
    P.zeros();
    for (int i = 0; i < 6; i++) {
        P += plaquetteStaples[i];
    }
    return (P.mat[0] + P.mat[8] + P.mat[16])/m_multiplicationFactor;
}

