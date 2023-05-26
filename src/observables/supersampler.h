/*!
 * \class SuperSampler
 *
 * \brief Observable for calculating the Weinberg operator, plaquette, energy and euclidean time topological charge.
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
 * See e.g. https://arxiv.org/abs/1711.04730 for a definition or the Weinberg operator.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef SUPERSAMPLER_H
#define SUPERSAMPLER_H

#include "math/lattice.h"
#include "observables/observables.h"
#include <map>

class SuperSampler : public Correlator
{
private:
    const std::string m_observableName = "Weinberg operator";
    const std::string m_observableNameCompact = "weinberg";

    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor, m_wMultiplicationFactor;
    double m_topCharge, m_energy, m_plaquette, m_weinberg;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_fieldTensorG[6];
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Container for the spatial observables
    std::vector<double> m_topctGatherVector;
    std::vector<double> m_wtGatherVector;
    std::vector<double> m_tempEucl;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
    ObservableStorer * m_topctObservable = nullptr;
    ObservableStorer * m_wObservable = nullptr; // Weinberg
    ObservableStorer * m_wtObservable = nullptr;

    std::map<int, std::map<int, int>> m_indexMap = {
        {0, {{1, 0}, {2, 1}, {3, 2}}},
        {1, {{2, 4}}},
        {3, {{1, 3}}},
        {2, {{3, 5}}},
    };

    inline int next_index(const int i)
    {
        /*Function for getting the next index.*/
        return i % 3 + 1;
    }

public:
    SuperSampler(const bool flow);
    ~SuperSampler();
    void calculate(Lattice<SU3> * lattice, const unsigned int iObs);
    void initializeObservableStorer(const bool storeFlowObservable);

    void writeObservableToFile(const double acceptanceRatio);
    void writeFlowObservablesToFile(const unsigned int configNumber);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(const unsigned int iObs);
    void printStatistics();
    std::string getObservableName() { return m_observableName; }
    std::vector<double> getObservablesVector(const unsigned int iObs);
    void copyObservable(const unsigned int iObs, const std::vector<double> &obs);
};

#endif // SUPERSAMPLER_H
