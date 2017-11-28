#ifndef TOPOLOGICALCHARGE_H
#define TOPOLOGICALCHARGE_H

#include "correlator.h"

class TopologicalCharge : public Correlator
{
private:
    const std::string m_observableName = "Topological Charge";
    const std::string m_observableNameCompact = "topc";
    double topCharge;
    double m_multiplicationFactor;
public:
    TopologicalCharge(bool storeFlowObservable);
    ~TopologicalCharge();

    void calculate(SU3 *clovers, int iObs);
    void calculate(Links *lattice, int iObs);

    // Stats
    void runStatistics();
    // Printers
    void printStatistics();
    // Getters
    std::string getObservableName() { return m_observableName; }
};

#endif // TOPOLOGICALCHARGE_H
