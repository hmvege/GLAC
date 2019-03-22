/*!
 * \class EnergyDensity
 *
 * \brief Solely calculates the energy density as given by,
 * \f[
 *      E = -\frac{1}{64|\Lambda|} \sum_{n\in\Lambda}\sum_{\mu<\nu} \big(\tilde{C}_{\mu\nu}(n)\big)^2.
 * \f]
 * The \f$F_{\mu\nu}\f$ is the clover.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef ENERGYDENSITY_H
#define ENERGYDENSITY_H

#include "correlator.h"
#include "math/functions.h"

class EnergyDensity : public Correlator
{
private:
    const std::string m_observableName = "Energy density";
    const std::string m_observableNameCompact = "energy";

    // Indexes for retrieving the clovers
    int mu = 0, rho, sigma;
    double m_energyDensity;
    double m_multiplicationFactor;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;
public:
    EnergyDensity(bool storeFlowObservable);
    ~EnergyDensity();
    void calculate(Lattice<SU3> *lattice, unsigned int iObs);

    // Printers
    void printStatistics();
    // Getters
    std::string getObservableName() { return m_observableName; }
    void runStatistics();
};

#endif // ENERGYDENSITY_H
