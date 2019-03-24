/*!
 * \class LatticeActionChargeDensity
 *
 * \brief Method for calculating and storing the topological charge and energy density at every point of the lattice.
 *
 * The energy is given by,
 * \f[
 *      E(n) = -\frac{1}{64} \sum_{\mu<\nu} \big(\tilde{C}_{\mu\nu}(n)\big)^2,
 * \f]
 * and the topological charge is given by,
 * \f[
 *      q(n) = \frac{1}{32\pi^2} \epsilon_{\mu\nu\rho\sigma} \mathrm{tr}\left[F_{\mu\nu}(n)F_{\rho\sigma}(n)\right].
 * \f]
 * where \f$n\f$ is the position in the lattice. \f$F_{\nu\mu}\f$ is the clover.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef LATTICEACTIONCHARGEDENSITY_H
#define LATTICEACTIONCHARGEDENSITY_H

#include "math/lattice.h"
#include "observables/observables.h"

class LatticeActionChargeDensity : public Correlator
{
private:
    int m_samplingFrequency = 25;
    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    double m_plaquette;
    Lattice <double> m_tempDiag, m_topCharge, m_energy;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
public:
    LatticeActionChargeDensity(bool flow);
    ~LatticeActionChargeDensity();
    void calculate(Lattice<SU3> * lattice, unsigned int iObs);
    void initializeObservableStorer(bool storeFlowObservable);

    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(unsigned int iFlow);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(unsigned int iObs);
    void printStatistics();
    std::vector<double> getObservablesVector(unsigned int iObs);
    void copyObservable(unsigned int iObs, std::vector<double> obs);
};

#endif // LATTICEACTIONCHARGEDENSITY_H
