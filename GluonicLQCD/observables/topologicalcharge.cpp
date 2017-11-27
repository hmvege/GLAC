#include "topologicalcharge.h"
#include "math/functions.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include "clover.h"
#include <cmath>

TopologicalCharge::TopologicalCharge(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_multiplicationFactor = 1.0/(16*16*M_PI*M_PI);
    populateLC(); // Fills the levi civita vector
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(false);
}

TopologicalCharge::~TopologicalCharge()
{

}

void TopologicalCharge::calculate(Links *lattice, int iObs)
{
    /*
     * Function to be used when no clover is provided. SHOULD BE TESTED
     */
    Clover Clov(m_storeFlowObservable);
    Clov.setN(m_N);
    Clov.setLatticeSize(m_latticeSize);
    topCharge = 0;

    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    Clov.calculateClover(lattice,i,j,k,l);
//                    for (unsigned int i = 0; i < m_leviCivita.size(); i++)
                    for (unsigned int i = 0; i < 3; i++)
                    {
                        traceRealMultiplication(Clov.m_clovers[2*i],Clov.m_clovers[2*i+1]);
//                        G1 = Clov.m_clovers[m_leviCivita[i].ci[0]];
//                        G2 = Clov.m_clovers[m_leviCivita[i].ci[1]];
////                        topCharge += traceSparseImagMultiplication(G1,G2)*m_leviCivita[i].sgn; // When off diagonal complex elements are zero
//                        topCharge += traceImagMultiplication(G1,G2)*m_leviCivita[i].sgn; // When off diagonal complex elements are not zero
                    }
                }
            }
        }
    }
//    return topCharge*m_multiplicationFactor;
    m_observable->m_observables[iObs] = topCharge*m_multiplicationFactor;
}

void TopologicalCharge::calculate(SU3 *clovers, int iObs)
{
    topCharge = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
//        G1 = clovers[m_leviCivita[i].ci[0]];
//        G2 = clovers[m_leviCivita[i].ci[1]];
//        topCharge += traceSparseImagMultiplication(G1,G2)*m_leviCivita[i].sgn;
//        topCharge += traceImagMultiplication(G1,G2)*m_leviCivita[i].sgn;
        topCharge += traceRealMultiplication(clovers[2*i],clovers[2*i+1]);
    }
    m_observable->m_observables[iObs] += topCharge*m_multiplicationFactor;
    // Jacks data should give for flow at t=0 3.29640613544198
//    Parallel::Communicator::setBarrier();
//    Parallel::Communicator::MPIExit("Exiting at wilson gauge action");

//    return topCharge*m_multiplicationFactor;
}

void TopologicalCharge::populateLC()
{
    int muNuOverCounter = 0;
    int rhoSigmaOverCounter = 0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            if (nu==mu) {
                muNuOverCounter++; // Acounts for overcounting
                continue;
            }

            for (int rho = 0; rho < 4; rho++) {
                if (rho==mu || rho==nu) {
                    rhoSigmaOverCounter++; // Acounts for overcounting
                    continue;
                }

                for (int sigma = 0; sigma < 4; sigma++) {

                    if (sigma==rho) {
                        rhoSigmaOverCounter++; // Acounts for overcounting
                        continue;
                    }
                    if (sigma==mu || sigma==nu) continue;
                    LeviCivita LC;
                    LC.setLC(mu, nu, rho, sigma);
                    LC.sgn = getLCSign(LC);
                    LC.setCI(cloverIndex(mu,nu - muNuOverCounter),cloverIndex(rho,sigma - rhoSigmaOverCounter));
                    m_leviCivita.push_back(LC);
                }
            }
            rhoSigmaOverCounter = 0;
        }
    }

    if (m_leviCivita.size() != 24)
    {
        cout << "Error: number of levi civita combinations is " << m_leviCivita.size() << ", not 24!" << endl;
        exit(1);
    }
}

int TopologicalCharge::getLCSign(LeviCivita LC)
{
    int sign = 1;
    int x = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            x = (LC.lc[i] - LC.lc[j]);
            sign *= (x > 0) - (x < 0);
        }
    }
    return sign;
}

void TopologicalCharge::printStatistics()
{
    m_observable->printStatistics();
}

void TopologicalCharge::runStatistics()
{
    /*
     * Statistics. Should perhaps make into its own class?
     */
    m_observable->gatherResults();
    m_observable->runStatistics();
//    int NObs = m_observable->m_NObs;
//    // Gathers results from processors
//    Parallel::Communicator::gatherDoubleResults(m_observable->m_observables,NObs);
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
