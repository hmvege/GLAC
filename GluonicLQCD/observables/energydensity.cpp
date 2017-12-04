#include "energydensity.h"
#include "clover.h"

EnergyDensity::EnergyDensity(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(false);
}

EnergyDensity::~EnergyDensity()
{

}

EnergyDensity::EnergyDensity(bool storeFlowObservable, double a, int latticeSize) : Correlator(storeFlowObservable)
{
    setLatticeSpacing(a);
    setLatticeSize(latticeSize);
//    m_multiplicationFactor = (m_a*m_a*m_a*m_a)/(3*double(m_latticeSize));
    m_multiplicationFactor = 1.0/double(m_latticeSize);
}

void EnergyDensity::setLatticeSpacing(double a) // NEED TO DOUBLE CHECK THIS WITH ANDREA!
{
//    m_multiplicationFactor = (a*a*a*a)/(3*double(m_latticeSize));
    m_a = a;
    m_multiplicationFactor = 1.0/double(m_latticeSize);
}

void EnergyDensity::setLatticeSize(int latticeSize)
{
    m_latticeSize = latticeSize;
    m_multiplicationFactor = 1.0/double(latticeSize);
}

void EnergyDensity::calculate(SU3 *clovers, int iObs)
{
    m_actionDensity = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
        m_actionDensity += traceRealMultiplication(clovers[2*i],clovers[2*i]); // Might check this one with Andrea
        m_actionDensity += traceRealMultiplication(clovers[2*i+1],clovers[2*i+1]); // Might check this one with Andrea
    }
    m_observable->m_observables[iObs] += m_actionDensity*m_multiplicationFactor; // Normalizes by volume
//    return m_actionDensity;//*m_multiplicationFactor; // Correct or not?
}

void EnergyDensity::calculate(Links *lattice, int iObs)
{
    // When clover is not provided
    Clover Clov(m_storeFlowObservable);
    Clov.setN(m_N);
    Clov.setLatticeSize(m_latticeSize);
    m_actionDensity = 0;
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
                        m_actionDensity += traceRealMultiplication(Clov.m_clovers[2*i],Clov.m_clovers[2*i]); // Might check this one with Andrea
                        m_actionDensity += traceRealMultiplication(Clov.m_clovers[2*i+1],Clov.m_clovers[2*i+1]); // Might check this one with Andrea
//                        m_actionDensity += traceSparseImagMultiplication(Clov.m_clovers[2*i],Clov.m_clovers[2*i+1]);
                    }
                }
            }
        }
    }
    m_observable->m_observables[iObs] = m_actionDensity*m_multiplicationFactor;
//    return m_actionDensity;//*m_multiplicationFactor; // Temporary off
}

void EnergyDensity::runStatistics()
{
    m_observable->gatherResults();
    m_observable->runStatistics();
}

void EnergyDensity::printStatistics()
{
    m_observable->printStatistics();
}
