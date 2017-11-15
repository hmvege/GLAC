#include "correlator.h"

Correlator::Correlator(bool storeFlowObservable)
{
    // Initiates the lattice dimensions
    m_N = new unsigned int[4];
    // Sets position vector to zero
    m_position = std::vector<int>(4,0);

    // Sets the lorentz indices to zero
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
    }
    storeFlow(storeFlowObservable);
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

double Correlator::getObservable(int iObs)
{
    return m_observable->getObservable(iObs);
}

void Correlator::writeStatisticsToFile(int iConfig)
{
//    printf("\nFunction for writing statistics to file not implemented for base correlator class!");
    m_observable->runStatistics();
    if (Parameters::getVerbose()) {
        m_observable->printStatistics();
    }
    m_observable->writeFlowObservableToFile(iConfig);
}

void Correlator::writeStatisticsToFile()
{
//    printf("\nFunction for writing statistics to file not implemented for base correlator class!");
    m_observable->runStatistics();
    if (Parameters::getVerbose()) {
        m_observable->printStatistics();
    }
    m_observable->writeObservableToFile();
}

void Correlator::storeFlow(bool storeFlowObservable)
{
    m_storeFlowObservable = storeFlowObservable;
    if (m_storeFlowObservable) {
        m_observable = new ObservableStorer(Parameters::getFlowSamplePoints());
    } else {
        m_observable = new ObservableStorer(Parameters::getConfigSamplePoints());
    }
}

void Correlator::printStatistics()
{
    if (Parameters::getVerbose()) m_observable->printStatistics();
}
