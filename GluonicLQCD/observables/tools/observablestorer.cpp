#include "observablestorer.h"
#include "io/observablesio.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include <cmath>
#include <mpi.h>

ObservableStorer::ObservableStorer(int NSize)
{
    m_NObs = NSize;
    // Initializes arrays
    m_observables = new double[m_NObs];
    m_observablesSquared = new double[m_NObs];
    for (int iObs = 0; iObs < m_NObs; iObs++) m_observables[iObs] = 0;
//    printf("\nSize off storeage array: %d. ProcessID: %d. ObsName: %s\n",m_NObs,Parallel::Communicator::getProcessRank(),m_observableName.c_str());
}

ObservableStorer::~ObservableStorer()
{
    delete [] m_observables;
    delete [] m_observablesSquared;
}

void ObservableStorer::gatherResults()
{
    // Performing an average over the Monte Carlo obtained values
    for (int iBuffer = 0; iBuffer < m_NObs; iBuffer++) {
        MPI_Allreduce(&m_observables[iBuffer], &m_observables[iBuffer], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    // Normalizing by the number of processors if specified
    if (m_normalizeObservableByProcessor) {
        double numprocs = double(Parallel::Communicator::getNumProc());
        for (int i = 0; i < m_NObs; i++) {
            m_observables[i] /= numprocs;
        }
    }
}

void ObservableStorer::runStatistics()
{
    // Normalizes by the number of processors
    double *observablesSquared = new double[m_NObs];
    if (m_normalizeObservableByProcessor) {
        for (int i = 0; i < m_NObs; i++) {
            observablesSquared[i] = m_observables[i]*m_observables[i];
        }
    }
    // Gets average of the observable
    for (int i = 0; i < m_NObs; i++)
    {
        m_averagedObservable += m_observables[i];
        m_averagedObservableSquared += observablesSquared[i];
    }
    m_averagedObservable /= double(m_NObs);
    m_averagedObservableSquared /= double(m_NObs);
    m_varianceObservable = (m_averagedObservableSquared - m_averagedObservable*m_averagedObservable)/double(m_NObs);
    m_stdObservable = sqrt(m_varianceObservable);
}

void ObservableStorer::printStatistics()
{
    if (Parallel::Communicator::getProcessRank() == 0)
    {
        printf("\n%-20s ", m_observableName.c_str());
        printf("%-20.15f ", m_averagedObservable);
        printf("%-20.15f ", m_varianceObservable);
        printf("%-20.15f ", m_stdObservable);
    }
}

void ObservableStorer::writeObservableToFile(double acceptanceRatio)
{
    IO::writeObservablesToFile(acceptanceRatio,m_averagedObservable,m_varianceObservable,m_stdObservable,m_observables,m_NObs,m_observableName);
}

void ObservableStorer::writeFlowObservableToFile(int configNumber)
{
    IO::writeFlowObservableToFile(m_observables,m_observableName,configNumber);
}
