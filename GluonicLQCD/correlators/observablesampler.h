#ifndef OBSERVABLESAMPLER_H
#define OBSERVABLESAMPLER_H

#include "correlator.h"
#include "plaquette.h"
#include "clover.h"
#include "energydensity.h"
#include "topologicalcharge.h"
#include "observablesampler.h"

class ObservableSampler : public Correlator
{
private:
    double m_a;
    double m_P, m_E, m_Q;
    Clover m_clover;
    Plaquette m_plaquette;
    EnergyDensity m_energyDensity;
    TopologicalCharge m_topologicalCharge;
public:
    ObservableSampler(double a, int latticeSize);
    ~ObservableSampler();

    double calculate(Links *lattice);
    void setLatticeSize(int latticeSize);
};

#endif // ObservableSampler_H
