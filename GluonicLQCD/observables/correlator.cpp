#include "correlator.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"

Correlator::Correlator(bool storeFlowObservable)
{
    storeFlow(storeFlowObservable);
    // Initiates the lattice dimensions
    m_a = Parameters::getLatticeSpacing();
    m_N = Parameters::getN();
    m_latticeSize = double(Parameters::getSubLatticeSize());
}

//void Correlator::setLatticeSize(int latticeSize) // MOVE TO CONSTRUCTOR?!
//{
//    m_latticeSize = double(latticeSize);
//}

Correlator::~Correlator()
{
    delete m_observable;
}

void Correlator::calculate(Lattice<SU3> *lattice, int iObs)
{
    /*
     * Default correlator is not implemented. Pushes to observable array at position iObs.
     */
    printf("\nENTERING BASE");
    lattice[0][iObs].zeros(); // TEMP
    printf("\nIf you see this, something is wrong! Should not call correlator.cpp");
}

//void Correlator::setN(std::vector<unsigned int> N) // MOVE INTO CONSTRUCTOR?
//{
//    for (int i = 0; i < 4; i++) {
//        m_N[i] = N[i];
//    }
//}

void Correlator::printHeader()
{
    printf("%-*s",m_headerWidth,m_observableName.c_str());
}

double Correlator::getObservable(int iObs)
{
    return (*m_observable)[iObs];
}

void Correlator::printObservable(int iObs)
{
    if (!m_storeFlowObservable) {
        printf("%-*.8f",m_headerWidth,(*m_observable)[iObs]);
    } else {
        printf("\n    %-*.8f",m_headerWidth,(*m_observable)[iObs]);
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
    m_observable->reset();
}

void Correlator::copyObservable(int iObs, std::vector<double> obs) {
    /*
     * Used when we already have calculated the observable in the
     */
    m_observable[iObs] = obs[0];
}

void Correlator::setObservable(int iObs, double obs) {
    m_observable[iObs] = obs;
}

std::vector<double> Correlator::getObservablesVector(int iObs) {
    std::vector<double> obs(1);
    obs[0] = (*m_observable)[iObs];
    return obs;
}

void Correlator::printStatistics()
{
    m_observable->printStatistics();
}
