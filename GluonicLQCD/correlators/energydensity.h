#ifndef ENERGYDENSITY_H
#define ENERGYDENSITY_H

#include "correlator.h"

class EnergyDensity : public Correlator
{
private:
    double m_multiplicationFactor;
public:
    SU3 m_clover[12];
    EnergyDensity(double a, double latticeVolume);
    void setClover(SU3 *clover);
    double calculate();
    double calculate(Links *lattice);
};

#endif // ENERGYDENSITY_H
