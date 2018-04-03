#ifndef TOPOLOGICALCHARGE_H
#define TOPOLOGICALCHARGE_H

#include "correlator.h"

class TopologicalCharge : public Correlator
{
private:
    const std::string m_observableName = "Topological Charge";
    const std::string m_observableNameCompact = "topc";

    // Indexes for retrieving the clovers
    int mu = 0, rho, sigma;
    double m_topCharge;
    double m_multiplicationFactor;

    Lattice <double> m_tempDiag;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;
public:
    TopologicalCharge(bool storeFlowObservable);
    ~TopologicalCharge();

    void calculate(Lattice<SU3> *lattice, unsigned int iObs);

    // Stats
    void runStatistics();
    // Printers
    void printStatistics();
    // Getters
    std::string getObservableName() { return m_observableName; }
};

#endif // TOPOLOGICALCHARGE_H
