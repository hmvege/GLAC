#include "plaquette.h"
#include "math/functions.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include <vector>
#include <cmath>

Plaquette::Plaquette(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(true);
    setLatticeSize(Parameters::getSubLatticeSize());
}

Plaquette::~Plaquette()
{
}

void Plaquette::setLatticeSize(unsigned long latticeSize)
{
    m_latticeSize = double(latticeSize);
    m_temp.allocate(m_N);
    m_multiplicationFactor = 1/(18.0*m_latticeSize); // 3 from SU3, 6 from number of plaquettes, 3*6=18
}

void Plaquette::calculate(Lattice<SU3> *lattice, unsigned int iObs)
{
    m_tempObservable = 0;
    for (unsigned int mu = 0; mu < 4; mu++) {
        for (unsigned int nu = mu+1; nu < 4; nu++) {
            m_temp = lattice[mu];
            m_temp *= shift(lattice[nu],FORWARDS,mu);
            m_temp *= inv(shift(lattice[mu],FORWARDS,nu));
            m_temp *= inv(lattice[nu]);
            m_tempObservable += sumRealTrace(m_temp);
        }
    }
    m_tempObservable *= m_multiplicationFactor;
    (*m_observable)[iObs] = m_tempObservable;
}

void Plaquette::runStatistics()
{
    /*
     * Statistics. Should perhaps make into its own class?
     */
    // Gathers results from processors
    m_observable->gatherResults();
    // Runs statistics, normalization defined in initialization
    m_observable->runStatistics();
}

void Plaquette::printObservable(unsigned int iObs)
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        if (!m_storeFlowObservable) {
            printf("%-*.8f", m_headerWidth, (*m_observable)[iObs]);
        } else {
            printf("\n%-4d %-*.8f", iObs, m_headerWidth, (*m_observable)[iObs]);
        }
    }
}

void Plaquette::printHeader()
{
    printf("%-*s", m_headerWidth, m_observableName.c_str());
}
