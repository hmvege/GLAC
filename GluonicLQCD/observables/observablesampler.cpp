#include "observablesampler.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"

ObservableSampler::ObservableSampler(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_clover = new Clover(storeFlowObservable);
    m_topologicalCharge = new TopologicalCharge(storeFlowObservable);
    m_energyDensity = new EnergyDensity(storeFlowObservable);
    m_plaquette = new Plaquette(storeFlowObservable);
    m_a = Parameters::getLatticeSpacing();
    m_energyDensity->setLatticeSpacing(m_a);
    setLatticeSize(Parameters::getSubLatticeSize());
    Parameters::getN(m_N);
    m_position = std::vector<int>(4,0);
}

ObservableSampler::~ObservableSampler()
{

}

void ObservableSampler::calculate(Links *lattice, int iObs)
{
    /*
     * Samples plaquette, topological charge and action/energy density.
     */
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    // Gets clover
                    m_clover->calculateClover(lattice,i,j,k,l);
                    m_plaquette->calculate(m_clover->m_plaquettes,iObs);
                    m_energyDensity->calculate(m_clover->m_clovers,iObs);
                    m_topologicalCharge->calculate(m_clover->m_clovers,iObs);
                }
            }
        }
    }
}

void ObservableSampler::setLatticeSize(int latticeSize)
{
    /*
     * Sets the lattice size of all the internal correlators.
     */
    m_latticeSize = latticeSize;
    m_clover->setLatticeSize(m_latticeSize);
    m_plaquette->setLatticeSize(m_latticeSize);
    m_energyDensity->setLatticeSize(m_latticeSize);
    m_topologicalCharge->setLatticeSize(m_latticeSize);
}

void ObservableSampler::writeStatisticsToFile(double acceptanceRatio)
{
    /*
     * Used by for writing the config stats to file.
     */
    m_plaquette->writeStatisticsToFile(acceptanceRatio);
    m_energyDensity->writeStatisticsToFile(acceptanceRatio);
    m_topologicalCharge->writeStatisticsToFile(acceptanceRatio);
}

void ObservableSampler::writeFlowObservablesToFile(int iFlow)
{
    /*
     * Used in flow.
     */
    m_plaquette->writeFlowObservablesToFile(iFlow);
    m_energyDensity->writeFlowObservablesToFile(iFlow);
    m_topologicalCharge->writeFlowObservablesToFile(iFlow);
}

void ObservableSampler::runStatistics()
{
    /*
     * Runs statistics for all the correlators.
     */
    m_plaquette->runStatistics();
    m_energyDensity->runStatistics();
    m_topologicalCharge->runStatistics();
}

void ObservableSampler::reset()
{
    /*
     * Resets internal values of the correlators(overflow otherwise).
     */
    m_plaquette->reset();
    m_energyDensity->reset();
    m_topologicalCharge->reset();
}

double ObservableSampler::getObservable(int iObs)
{
    /*
     * Returns plaquette value only.
     */
    return m_observable->getObservable(iObs);
}

void ObservableSampler::printHeader()
{
    printf("%-*s %-*s %-*s",
           m_headerWidth,m_plaquette->getObservableName().c_str(),
           m_headerWidth,m_topologicalCharge->getObservableName().c_str(),
           m_headerWidth,m_energyDensity->getObservableName().c_str());
}

void ObservableSampler::printObservable(int iObs)
{
    printf("%-*.4f %-*.4f %-*.4f",
           m_headerWidth,m_plaquette->getObservable(iObs),
           m_headerWidth,m_topologicalCharge->getObservable(iObs),
           m_headerWidth,m_energyDensity->getObservable(iObs));
}
