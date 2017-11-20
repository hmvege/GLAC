#ifndef ENERGYDENSITY_H
#define ENERGYDENSITY_H

#include "correlator.h"

class EnergyDensity : public Correlator
{
private:
    double m_multiplicationFactor;
    double m_actionDensity;
    const std::string m_observableName = "Energy density";
    const std::string m_observableNameCompact = "energy";
public:
    EnergyDensity(bool storeFlowObservable, double a, int latticeSize);
    EnergyDensity(bool storeFlowObservable);
    ~EnergyDensity();
    void calculate(SU3 *clovers, int iObs);
    void calculate(Links *lattice, int iObs);

    // Printers
    void printStatistics();
    // Setters
    void setLatticeSpacing(double a);
    // Getters
    virtual std::string getObservableName() { return m_observableName; }
};

#endif // ENERGYDENSITY_H
