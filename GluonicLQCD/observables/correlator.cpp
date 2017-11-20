#include "correlator.h"

Correlator::Correlator(bool storeFlowObservable)
{
    storeFlow(storeFlowObservable);
    // Initiates the lattice dimensions
    m_N = new unsigned int[4];
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
    delete [] m_N;
}

void Correlator::calculate(Links * lattice, int iObs)
{
    /*
     * Default correlator is not implemented. Pushes to observable array at position iObs.
     */
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

void Correlator::setN(unsigned int *N) // MOVE INTO CONSTRUCTOR?
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
    printf("\n%-4d %-*.8f",iObs,m_headerWidth,m_observable->m_observables[iObs]);
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
        m_observable = new ObservableStorer(Parameters::getNFlows());
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
