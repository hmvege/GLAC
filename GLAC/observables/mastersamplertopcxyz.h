/*!
 * \class MasterSamplerTopcXYZ
 *
 * \brief Observable for calculating the plaquette, energy and euclidean time topological charge.
 *
 * Plaqeutte definition:
 * \f[
 *      P = \frac{1}{16|\Lambda|} \sum_{n\in\Lambda} \sum_{\mu < \nu} \Re \mathrm{tr} P_{\mu\nu}.
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
 *      q(n_e) = \frac{1}{32\pi^2} \sum_{n\in N^3} \epsilon_{\mu\nu\rho\sigma} \mathrm{tr}\left[F_{\mu\nu}(n)F_{\rho\sigma}(n)\right],
 * \f]
 * where \f$n_e\f$ is the position in the lattice in Euclidean time. \f$F_{\nu\mu}\f$ is the clover.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef MASTERSAMPLERTOPCXYZ_H
#define MASTERSAMPLERTOPCXYZ_H

#include "math/lattice.h"
#include "observables/observables.h"

class MasterSamplerTopcXYZ : public Correlator
{
private:
    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    double m_topCharge, m_energy, m_plaquette;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Container for the topc xyz observable
    std::vector<double> m_topctGatherVector;//m_tempTopctArray;
    std::vector<double> m_tempTopct;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
    ObservableStorer * m_topctObservable = nullptr;
public:
    MasterSamplerTopcXYZ(bool flow);
    ~MasterSamplerTopcXYZ();
    void calculate(Lattice<SU3> * lattice, unsigned int iObs);
    void initializeObservableStorer(bool storeFlowObservable);

    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(unsigned int configNumber);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(unsigned int iObs);
    void printStatistics();
    std::vector<double> getObservablesVector(unsigned int iObs);
    void copyObservable(unsigned int iObs, std::vector<double> obs);
};

#endif // MASTERSAMPLERTOPCXYZ_H
