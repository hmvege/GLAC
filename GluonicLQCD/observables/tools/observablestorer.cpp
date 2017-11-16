#include "observablestorer.h"

ObservableStorer::ObservableStorer(int NSize)
{
    // Initializes arrays
    m_NObs = NSize;
    m_observables = new double[m_NObs];
    for (int iObs = 0; iObs < m_NObs; iObs++) m_observables[iObs] = 0;
//    printf("\nSize off storeage array: %d. ProcessID: %d. ObsName: %s\n",m_NObs,Parallel::Communicator::getProcessRank(),m_observableName.c_str());
}

ObservableStorer::~ObservableStorer()
{
    delete [] m_observables;
}


void ObservableStorer::runStatistics()
{
    // Performing an average over the Monte Carlo obtained values
    for (int iBuffer = 0; iBuffer < m_NObs; iBuffer++) {
        MPI_Allreduce(&m_observables[iBuffer], &m_observables[iBuffer], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    // Normalizes by the number of processors
    double *observablesSquared = new double[m_NObs];
//    double numprocs = double(Parallel::Communicator::getNumProc());
//    if (m_normalizeObservableByProcessor) {
//        for (int i = 0; i < m_NObs; i++) {
//            m_observables[i] /= numprocs;
//            observablesSquared[i] = m_observables[i]*m_observables[i];
//        }
//    }
    // Gets average of the observable
    double averagedObservableSquared = 0;
    for (int i = 0; i < m_NObs; i++)
    {
        m_averagedObservable += m_observables[i];
        averagedObservableSquared += observablesSquared[i];
    }
    m_averagedObservable /= double(m_NObs);
    averagedObservableSquared /= double(m_NObs);
    m_varianceObservable = (averagedObservableSquared - m_averagedObservable*m_averagedObservable)/double(m_NObs);
    m_stdObservable = sqrt(m_varianceObservable);
    delete [] observablesSquared;
}

void ObservableStorer::printStatistics()
{
    if (Parallel::Communicator::getProcessRank() == 0)
    {
//        printf("\nObservable: %s", m_observableName.c_str());
//        printf("\nAverage: %-20.15f", m_averagedObservable);
//        printf("\nVariance%-20.15f", m_varianceObservable);
//        printf("\nStandard deviation: %-20.15f", m_stdObservable);
        printf("\n%-20s ", m_observableName.c_str());
        printf("%-20.15f", m_averagedObservable);
        printf("%-20.15f", m_varianceObservable);
        printf("%-20.15f", m_stdObservable);
    }
}

void ObservableStorer::writeObservableToFile()
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        IO::writeDataToFile(m_averagedObservable,m_varianceObservable,m_stdObservable,m_observables,m_observableName);
    }
}

void ObservableStorer::writeFlowObservableToFile(int configNumber)
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        IO::writeFlowObservableToFile(m_averagedObservable,m_varianceObservable,m_stdObservable,m_observables,m_observableName,configNumber);
    }
}
