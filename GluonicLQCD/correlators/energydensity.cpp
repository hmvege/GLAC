#include "energydensity.h"

EnergyDensity::EnergyDensity(double a, double latticeVolume) : Correlator()
{
    m_multiplicationFactor = (a*a*a*a)/(3*latticeVolume);
}

void EnergyDensity::setClover(SU3 *clover)
{
    for (int i = 0; i < 12; i++) {
        m_clover[i] = clover[i];
    }
}

double EnergyDensity::calculate()
{
    // When clover is provided
    return 0;
}

double EnergyDensity::calculate(Links *lattice)
{
    // When clover is not provided
    return 0;
}
