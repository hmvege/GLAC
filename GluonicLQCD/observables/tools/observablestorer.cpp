#include "observablestorer.h"
#include "io/observablesio.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include <cmath>
#include <mpi.h>

ObservableStorer::ObservableStorer(int NSize)
{
    /*
     * Initialization of a observable storage container.
     * Args:
     *  NSize: the number of observables we will sample.
     */
    m_NObs = NSize;
    // Initializes arrays for storing the observables
    m_observables = new double[m_NObs];
    m_observablesSquared = new double[m_NObs];
    for (int iObs = 0; iObs < m_NObs; iObs++) m_observables[iObs] = 0;
}

ObservableStorer::~ObservableStorer()
{
    delete [] m_observables;
    delete [] m_observablesSquared;
}

void ObservableStorer::gatherResults()
{
    /*
     * Gather all observable data from all of the processors
     */
    // Temporary buffer for summing the observables
    double tempBuffer[m_NObs];

    for (int iObs = 0; iObs < m_NObs; iObs++) {
        tempBuffer[iObs] = 0;
    }

    // Performing an average over the Monte Carlo obtained values
    MPI_Reduce(m_observables, tempBuffer, m_NObs, MPI_DOUBLE, MPI_SUM, 0, Parallel::ParallelParameters::ACTIVE_COMM);

    // Retrieving observable from temporary buffer
    for (int iBuffer = 0; iBuffer < m_NObs; iBuffer++) {
        m_observables[iBuffer] = tempBuffer[iBuffer];
    }

    // Normalizing by the number of processors if specified
    if (m_normalizeObservableByProcessor) {
        double numprocs = double(Parallel::Communicator::getNumProc());
        for (int i = 0; i < m_NObs; i++) {
            m_observables[i] /= double(numprocs);
        }
    }
}

void ObservableStorer::runStatistics()
{
    /*
     * Performs the statistics on the observables.
     * Best used in conjecture with ObservableStorer::gatherResults().
     */
    for (int i = 0; i < m_NObs; i++) {
        m_observablesSquared[i] = m_observables[i]*m_observables[i];
    }
    // Gets average of the observable
    for (int i = 0; i < m_NObs; i++)
    {
        m_averagedObservable += m_observables[i];
        m_averagedObservableSquared += m_observablesSquared[i];
    }
    m_averagedObservable /= double(m_NObs);
    m_averagedObservableSquared /= double(m_NObs);
    m_varianceObservable = (m_averagedObservableSquared - m_averagedObservable*m_averagedObservable)/double(m_NObs);
    m_stdObservable = sqrt(m_varianceObservable);
}

void ObservableStorer::printStatistics()
{
    /*
     * Prints the statistics from one processor.
     */
    if (Parallel::Communicator::getProcessRank() == 0)
    {
        printf("\n%-20s ", m_observableName.c_str());
        printf("%-25.15f ", m_averagedObservable);
        printf("%-25.15f ", m_varianceObservable);
        printf("%-25.15f ", m_stdObservable);
    }
}

void ObservableStorer::writeObservableToFile(double acceptanceRatio)
{
    /*
     * Writes a regular observable to file
     */
    IO::writeObservablesToFile(acceptanceRatio, m_averagedObservable, m_varianceObservable,
                               m_stdObservable, m_observables, m_NObs, m_observableName);
}

void ObservableStorer::writeFlowObservableToFile(int configNumber)
{
    /*
     * Writes a flow observable to file
     */
    IO::writeFlowObservableToFile(m_observables,m_observableName,configNumber);
}

void ObservableStorer::reset()
{
    /*
     * Sets all the observables to zero. A must-have when flowing.
     */
    for (int i = 0; i < m_NObs; i++)
    {
        m_observables[i] = 0;
        m_observablesSquared[i] = 0;
    }
    m_averagedObservable = 0;
    m_averagedObservableSquared = 0;
    m_varianceObservable = 0;
    m_stdObservable = 0;
}
