/*!
 * \class MasterSampler
 *
 * \brief Observable for calculating the plaquette, topological charge and energy.
 *
 * Plaqeutte definition:
 * \f[
 *      P = \frac{1}{16|\Lambda|} \sum_{n\in\Lambda} \sum_{\mu < \nu} \Re \mathrm{tr} P_{\mu\nu},
 * \f]
 * where
 * \f[
 *      P_{\mu\nu}=U_\mu(n) U_{\nu}(n+\hat{\mu}) U_{\mu}(n+\hat{\nu})^\dagger U_{\nu} (n)^\dagger
 * \f]
 * is the plquette.
 * Energy definition:
 * \f[
 *      E = -\frac{1}{64|\Lambda|} \sum_{n\in\Lambda} \sum_{\mu<\nu} \left(\tilde{C}_{\mu\nu}(n)\right)^2.
 * \f]
 * Topological charge definition:
 * \f[
 *      Q = \frac{1}{32\pi^2} \sum_{n\in\Lambda} \epsilon_{\mu\nu\rho\sigma} \mathrm{tr}\left[F_{\mu\nu}(n)F_{\rho\sigma}(n)\right],
 * \f]
 * where \f$n\f$ is the position in the lattice. \f$F_{\nu\mu}\f$ is the clover.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef MASTERSAMPLER_H
#define MASTERSAMPLER_H

#include "math/lattice.h"
#include "observables/observables.h"

class MasterSampler : public Correlator
{
private:
    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    double m_topCharge, m_energy, m_plaquette;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
public:
    MasterSampler(const bool flow);
    ~MasterSampler();
    void calculate(Lattice<SU3> * lattice, const unsigned int iObs);
    void initializeObservableStorer(const bool storeFlowObservable);

    void writeObservableToFile(const double acceptanceRatio);
    void writeFlowObservablesToFile(const unsigned int iFlow);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(const unsigned int iObs);
    void printStatistics();
    std::vector<double> getObservablesVector(const unsigned int iObs);
    void copyObservable(const unsigned int iObs, const std::vector<double> &obs);
};

#endif // MASTERSAMPLER_H
