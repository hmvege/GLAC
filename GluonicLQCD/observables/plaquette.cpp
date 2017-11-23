#include "plaquette.h"
#include "math/links.h"
#include "math/functions.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include <vector>
#include <cmath>

Plaquette::Plaquette(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(true);
}

Plaquette::~Plaquette()
{
}

void Plaquette::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
    m_multiplicationFactor = 18.0*m_latticeSize; // 3 from SU3, 6 from number of plaquettes, 3*6=18
}

void Plaquette::calculate(Links *lattice, int iObs)
{
    P.zeros();
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    for (int mu = 0; mu < 4; mu++) {
                        updateMuIndex(mu); // Inline function
                        for (int nu = mu+1; nu < 4; nu++) {
                            updateNuIndex(nu); // Inline function
                            PTemp = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
                            PTemp *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
                            PTemp *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
                            PTemp *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
                            P += PTemp;
                        }
                    }
                }
            }
        }
    }
    m_observable->m_observables[iObs] = (P.mat[0] + P.mat[8] + P.mat[16])/m_multiplicationFactor;
}

void Plaquette::calculate(SU3 *plaquetteStaples, int iObs)
{
    P.zeros();
    for (int i = 0; i < 6; i++) {
        P += plaquetteStaples[i];
    }
    m_observable->m_observables[iObs] += (P.mat[0] + P.mat[8] + P.mat[16])/m_multiplicationFactor;
}

void Plaquette::runStatistics()
{
    /*
     * Statistics. Should perhaps make into its own class?
     */
//    int NObs = m_observable->m_NObs;
    // Gathers results from processors
    m_observable->gatherResults();
    m_observable->runStatistics();
//    Parallel::Communicator::gatherDoubleResults(m_observable->m_observables,NObs);
//    for (int iObs = 0; iObs < NObs; iObs++) {
//        m_observable->m_observables[iObs] /= double(Parallel::Communicator::getNumProc()); // For plaquette only
//    }
//    // Temp holders
//    double averagedObservableSquared = 0;
//    for (int iObs = 0; iObs < NObs; iObs++) {
//        m_observable->m_averagedObservable += m_observable->m_observables[iObs];
//        averagedObservableSquared += m_observable->m_observables[iObs]*m_observable->m_observables[iObs];
//    }
//    averagedObservableSquared /= double(NObs);
//    m_observable->m_averagedObservable /= double(NObs);
//    m_observable->m_varianceObservable = (averagedObservableSquared - m_observable->m_averagedObservable*m_observable->m_averagedObservable)/double(NObs);
//    m_observable->m_stdObservable = sqrt(m_observable->m_varianceObservable);
}

void Plaquette::printObservable(int iObs)
{
    if (!m_storeFlowObservable) {
        printf("%-*.8f",m_headerWidth,m_observable->m_observables[iObs]);
    } else {
        printf("\n    %-*.8f",m_headerWidth,m_observable->m_observables[iObs]);
    }
}

void Plaquette::printHeader()
{
    printf("%-*s",m_headerWidth,m_observableName.c_str());
}
