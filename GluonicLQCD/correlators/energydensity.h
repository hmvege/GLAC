#ifndef ENERGYDENSITY_H
#define ENERGYDENSITY_H

#include "correlator.h"

class EnergyDensity : public Correlator
{
private:
    double m_multiplicationFactor;
    double m_actionDensity;
public:
    SU3 m_clover[12];
    SU3 m_temp;
    EnergyDensity(double a, int latticeSize);
    void setClover(SU3 *clover);
    double calculate();
    double calculate(Links *lattice);
};

#endif // ENERGYDENSITY_H
