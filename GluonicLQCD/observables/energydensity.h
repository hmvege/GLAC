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
