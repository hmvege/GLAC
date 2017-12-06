#include "plaquette.h"
#include "math/links.h"
#include "math/functions.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include <vector>
#include <cmath>

//using LatOps::FORWARDS;
//using LatOps::BACKWARDS;
//using LatOps::shift;
//using LatOps::sumRealTrace;
//using LatOps::inv;

Plaquette::Plaquette(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(true);
}

Plaquette::~Plaquette()
{
}

void Plaquette::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
    m_temp.allocate(m_N);
    m_multiplicationFactor = 1/(18.0*m_latticeSize); // 3 from SU3, 6 from number of plaquettes, 3*6=18
}

void Plaquette::calculate(Lattice<SU3> *lattice, int iObs)
{
    m_tempObservable = 0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu+1; nu < 4; nu++) {
            m_temp = lattice[mu];
            m_temp *= shift(lattice[nu],FORWARDS,mu);
            m_temp *= inv(shift(lattice[mu],FORWARDS,nu));
            m_temp *= inv(lattice[nu]);
            m_tempObservable += sumRealTrace(m_temp);
        }
    }
    m_tempObservable *= m_multiplicationFactor;
    m_observable[iObs] = m_tempObservable;
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

void Plaquette::printObservable(int iObs)
{
    if (!m_storeFlowObservable) {
        printf("%-*.8f",m_headerWidth,(*m_observable)[iObs]);
    } else {
        printf("\n    %-*.8f",m_headerWidth,(*m_observable)[iObs]);
    }
}

void Plaquette::printHeader()
{
    printf("%-*s",m_headerWidth,m_observableName.c_str());
}
