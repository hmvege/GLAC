#include "correlator.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"

Correlator::Correlator(bool storeFlowObservable)
{
    storeFlow(storeFlowObservable);
    // Initiates the lattice dimensions
    // Sets position vector to zero
    m_position = std::vector<int>(4,0);

    // Sets the lorentz indices to zero
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
}

void Correlator::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
}

Correlator::~Correlator()
{
    delete m_observable;
}

void Correlator::calculate(Links * lattice, int iObs)
{
    /*
     * Default correlator is not implemented. Pushes to observable array at position iObs.
     */
    printf("\nENTERING BASE");
    lattice[iObs].U[0].zeros(); // TEMP
    printf("\nIf you see this, something is wrong! Should not call correlator.cpp");
}

void Correlator::calculate(SU3 *U, int iObs)
{
    /*
     * Default correlator is not implemented when only given a SU3 matrix. Pushes to observable array at position iObs.
     */
    U[iObs % 4].zeros(); // TEMP
    printf("\nIf you see this, something is wrong! Should not call correlator.cpp");
}

void Correlator::setN(std::vector<unsigned int> N) // MOVE INTO CONSTRUCTOR?
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}

void Correlator::printHeader()
{
    printf("%-*s",m_headerWidth,m_observableName.c_str());
}

double Correlator::getObservable(int iObs)
{
    return m_observable->m_observables[iObs];
}

void Correlator::printObservable(int iObs)
{
    if (!m_storeFlowObservable) {
        printf("%-*.8f",m_headerWidth,m_observable->m_observables[iObs]);
    } else {
        printf("\n    %-*.8f",m_headerWidth,m_observable->m_observables[iObs]);
    }
}

void Correlator::runStatistics()
{
    m_observable->gatherResults();
    m_observable->runStatistics();
}

void Correlator::writeFlowObservablesToFile(int iFlow)
{
    m_observable->gatherResults();
    m_observable->writeFlowObservableToFile(iFlow);
}

void Correlator::writeStatisticsToFile(double acceptanceRatio)
{
    m_observable->writeObservableToFile(acceptanceRatio);
}

void Correlator::storeFlow(bool storeFlowObservable)
{
    m_storeFlowObservable = storeFlowObservable;
    if (m_storeFlowObservable) {
        m_observable = new ObservableStorer(Parameters::getNFlows() + 1); // +1 as we are storing the initial value at t=0 as well.
    } else {
        if (Parameters::getStoreThermalizationObservables()) {
            m_observable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
        } else {
            m_observable = new ObservableStorer(Parameters::getNCf());
        }
    }
}

void Correlator::reset()
{
    for (int i = 0; i < m_observable->m_NObs; i++) {
        m_observable->m_observables[i] = 0;
        m_observable->m_observablesSquared[i] = 0;
    }
    m_observable->m_stdObservable = 0;
    m_observable->m_varianceObservable = 0;
    m_observable->m_averagedObservable = 0;
    m_observable->m_averagedObservableSquared = 0;
}

void Correlator::copyObservable(int iObs, std::vector<double> obs) {
    /*
     * Used when we already have calculated the observable in the
     */
    m_observable->m_observables[iObs] = obs[0];
}

void Correlator::setObservable(int iObs, double obs) {
    m_observable->m_observables[iObs] = obs;
}

std::vector<double> Correlator::getObservablesVector(int iObs) {
    std::vector<double> obs(1);
    obs[0] = m_observable->getObservable(iObs);
    return obs;
}

void Correlator::printStatistics()
{
    m_observable->printStatistics();
}
