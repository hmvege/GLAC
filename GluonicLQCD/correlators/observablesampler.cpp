#include "observablesampler.h"

#include "parallelization/indexorganiser.h"

#include <iostream>
using std::cout;
using std::endl;

ObservableSampler::ObservableSampler(unsigned int *N, int latticeSize, double a, IndexOrganiser *Index)// : Correlator()
{
    m_a = a;
    setLatticeSize(latticeSize);
    m_energyDensity.setLatticeSpacing(m_a);
    setN(N);
    m_position = std::vector<int>(4,0);
    m_clover.initializeIndexHandler(Index);
    m_plaquette.initializeIndexHandler(Index);
    m_energyDensity.initializeIndexHandler(Index);
    m_topologicalCharge.initializeIndexHandler(Index);
}

ObservableSampler::~ObservableSampler()
{

}

void ObservableSampler::calculate(Links *lattice)
{
    /*
     * Samples plaquette, topological charge and action/energy density.
     */
    m_P = 0;
    m_E = 0;
    m_Q = 0;
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    // Gets clover
                    m_clover.calculateClover(lattice,i,j,k,l);
                    m_P += m_plaquette.calculate(m_clover.m_plaquettes);
                    m_E += m_energyDensity.calculate(m_clover.m_clovers);
                    m_Q += m_topologicalCharge.calculate(m_clover.m_clovers);
                }
            }
        }
    }
}

void ObservableSampler::setLatticeSize(int latticeSize)
{
    m_latticeSize = latticeSize;
    m_clover.setLatticeSize(m_latticeSize);
    m_plaquette.setLatticeSize(m_latticeSize);
    m_energyDensity.setLatticeSize(m_latticeSize);
    m_topologicalCharge.setLatticeSize(m_latticeSize);
}

void ObservableSampler::setN(unsigned int *N)
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}
