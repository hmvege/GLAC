#include "correlator.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"

/*!
 * \brief Correlator::Correlator
 * \param storeFlowObservable specifies if one is to sample (and store) a flow observable when initialized.
 */
Correlator::Correlator(bool storeFlowObservable)
{
    initializeObservableStorer(storeFlowObservable);
    // Initiates the lattice dimensions
    m_a = Parameters::getLatticeSpacing();
    m_N = Parameters::getN();
    m_latticeSize = double(Parameters::getSubLatticeSize());
}

/*!
 * \brief Correlator::Correlator
 *
 * When no storeFlowObservable is passed, default is to store observables for the configurations.
 */
Correlator::Correlator()
{
    // Initiates the lattice dimensions
    m_a = Parameters::getLatticeSpacing();
    m_N = Parameters::getN();
    m_latticeSize = double(Parameters::getSubLatticeSize());
}

Correlator::~Correlator()
{
    // Freeing observable storage container
    delete m_observable;
}

/*!
 * \brief Correlator::calculate calculates a observable
 * \param lattice pointer to lattice of SU3 matrices
 * \param iObs number of the observable
 */
void Correlator::calculate(Lattice<SU3> *lattice, unsigned int iObs)
{
    /*
     * Default correlator is not implemented. Pushes to observable array at position iObs.
     */
    lattice[0][iObs].zeros();
    printf("\nIf you see this, something is wrong! Should not call correlator.cpp");
}

/*!
 * \brief Correlator::printHeader prints a the observable name for output in header.
 */
void Correlator::printHeader()
{
    printf("%-*s",m_headerWidth,m_observableName.c_str());
}

/*!
 * \brief Correlator::getObservable
 * \param iObs number of the observable to print.
 * \return the observable for given iObs.
 */
double Correlator::getObservable(unsigned int iObs)
{
    return (*m_observable)[iObs];
}

/*!
 * \brief Correlator::printObservable prints the observable for the output.
 * \param iObs
 */
void Correlator::printObservable(unsigned int iObs)
{
    if (Parallel::Communicator::getProcessRank() == 0)
    {
        if (!m_storeFlowObservable)
        {
            printf("%-*.8f",m_headerWidth,(*m_observable)[iObs]);
        }
        else
        {
            printf("\n%-4d %-*.8f",iObs,m_headerWidth,(*m_observable)[iObs]);
        }
    }
}

/*!
 * \brief Correlator::runStatistics
 *
 * Gathers results and performs statistics in the ObservableStorer object.
 */
void Correlator::runStatistics()
{
    /*
     * Used before writeObservableToFile()
     */
    m_observable->gatherResults();
    m_observable->runStatistics();
}

/*!
 * \brief Correlator::writeFlowObservablesToFile
 * \param iFlow
 *
 * Gathers results and performs statistics in the ObservableStorer object.
 */
void Correlator::writeFlowObservablesToFile(unsigned int iFlow)
{
    m_observable->gatherResults();
    m_observable->writeFlowObservableToFile(iFlow);
}

/*!
 * \brief Correlator::writeObservableToFile calls the ObservableStorer for writing observable to file.
 * \param acceptanceRatio
 */
void Correlator::writeObservableToFile(double acceptanceRatio)
{
    m_observable->writeObservableToFile(acceptanceRatio);
}

/*!
 * \brief Correlator::initializeObservableStorer initializes ObservableStorer.
 * \param storeFlowObservable
 */
void Correlator::initializeObservableStorer(bool storeFlowObservable)
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

/*!
 * \brief Correlator::reset resets ObservableStorer.
 */
void Correlator::reset()
{
    m_observable->reset();
}

/*!
 * \brief Correlator::copyObservable copies the observable if has already been calculated at zero flow time.
 * \param iObs observable number.
 * \param obs vector containing observable information. Is a vector in case it consists of multiple values(e.g. topc in Euclidean time).
 */
void Correlator::copyObservable(unsigned int iObs, std::vector<double> obs) {
    /*
     * Used when we already have calculated the observable in the
     */
    (*m_observable)[iObs] = obs[0];
}

/*!
 * \brief Correlator::getObservablesVector
 * \param iObs
 * \return a vector containing the observables at given iObs.
 */
std::vector<double> Correlator::getObservablesVector(unsigned int iObs) {
    std::vector<double> obs(1);
    obs[0] = (*m_observable)[iObs];
    return obs;
}

/*!
 * \brief Correlator::printStatistics prints statistics from ObservableStorer.
 */
void Correlator::printStatistics()
{
    m_observable->printStatistics();
}
