#include "topologicalcharge.h"
#include "math/functions.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include "clover.h"
#include <cmath>

TopologicalCharge::TopologicalCharge(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_multiplicationFactor = 1.0/(16*16*M_PI*M_PI);
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(false);
}

TopologicalCharge::~TopologicalCharge()
{

}

void TopologicalCharge::calculate(Lattice<SU3> *lattice, int iObs)
{
    /*
     * Function to be used when no clover is provided. SHOULD BE TESTED
     */
    Clover Clov(m_storeFlowObservable);
    Clov.setN(m_N);
    Clov.setLatticeSize(m_latticeSize);
    topCharge = 0;

    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    Clov.calculateClover(lattice,i,j,k,l);
                    for (unsigned int i = 0; i < 3; i++)
                    {
                        topCharge -= traceRealMultiplication(Clov.m_clovers[2*i],Clov.m_clovers[2*i+1]);
                    }
                }
            }
        }
    }
    (*m_observable)[iObs] = topCharge*m_multiplicationFactor;
}

void TopologicalCharge::printStatistics()
{
    m_observable->printStatistics();
}

void TopologicalCharge::runStatistics()
{
    /*
     * Statistics. Should perhaps make into its own class?
     */
    m_observable->gatherResults();
    m_observable->runStatistics();
}
