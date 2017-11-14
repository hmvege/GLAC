#ifndef ENERGYDENSITY_H
#define ENERGYDENSITY_H

#include "correlator.h"

class EnergyDensity : public Correlator
{
private:
    double m_multiplicationFactor;
    double m_actionDensity;
public:
    EnergyDensity(double a, int latticeSize);
    EnergyDensity();
    ~EnergyDensity();
    void calculate(SU3 *clovers, int i);
    void calculate(Links *lattice, int i);
    void setLatticeSpacing(double a);
};

#endif // ENERGYDENSITY_H
