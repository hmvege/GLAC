#include "observablestorer.h"

ObservableStorer::ObservableStorer(int NSize, std::string observableName, bool procNormalize) :
    m_NSize(NSize),
    m_observableName(observableName),
    m_procNormalize(procNormalize)
{
    m_observables = new double[m_NSize];
    m_observablesSquared = new double[m_NSize];
}

ObservableStorer::~ObservableStorer()
{
    delete [] m_observables;
    delete [] m_observablesSquared;
}

void ObservableStorer::pushObservable(double newObs, int position)
{
    m_observables[position] = newObs;
    m_observablesSquared[position] = newObs*newObs;
}

void ObservableStorer::runStatistics()
{
    // Performing an average over the Monte Carlo obtained values
    MPI_Allreduce(&m_observables, &m_observables, m_NSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&m_observablesSquared, &m_observablesSquared, m_NSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double numprocs = double(Parallel::Communicator::getNumProc());
    double averagedObservableSquared = 0;
    // Normalizes by the number of processors
    if (m_procNormalize) {
        for (int i = 0; i < m_NSize; i++) {
            m_observables[i] /= numprocs;
            m_observablesSquared[i] /= numprocs;
        }
    }
    // Gets average of the observable
    for (int i = 0; i < m_NSize; i++)
    {
        m_averagedObservable += m_observables[i];
        averagedObservableSquared += m_observablesSquared[i];
    }
    m_averagedObservable /= double(m_NSize);
    averagedObservableSquared /= double(m_NSize);
    m_varianceObservable = (averagedObservableSquared - m_averagedObservable*m_averagedObservable)/double(m_NSize);
    m_stdObservable = sqrt(m_varianceObservable);
    if (Parallel::Communicator::getProcessRank() == 0 && Parameters::getVerbose())
    {
        printf("\nObservable: %s", m_observableName.c_str());
        printf("\nAverage: %-20.15f", m_averagedObservable);
        printf("\nStandard deviation: %-20.15f", m_averagedObservable);
        printf("\nVariance%-20.15f", m_averagedObservable);
    }
}

void ObservableStorer::writeObservableToFile()
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        IO::writeDataToFile(m_averagedObservable,m_varianceObservable,m_stdObservable,m_observables,m_observableName);
    }
}

void ObservableStorer::writeFlowObservableToFile(double * flowTime, int configNumber)
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        IO::writeFlowObservableToFile(m_averagedObservable,m_varianceObservable,m_stdObservable,flowTime,m_observables,m_observableName,configNumber);
    }
}
