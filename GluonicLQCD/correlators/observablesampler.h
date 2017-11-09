#ifndef OBSERVABLESAMPLER_H
#define OBSERVABLESAMPLER_H

#include "plaquette.h"
#include "clover.h"
#include "energydensity.h"
#include "topologicalcharge.h"
#include "observablesampler.h"

class ObservableSampler// : public Correlator
{
private:
    unsigned int m_N[4];
    int m_latticeSize;
    double m_a;
    double m_P, m_E, m_Q;
    std::vector<int> m_position;

    // All observables in program
    Clover m_clover;
    Plaquette m_plaquette;
    EnergyDensity m_energyDensity;
    TopologicalCharge m_topologicalCharge;
public:
    ObservableSampler(unsigned int *N, int latticeSize, double a);
    ~ObservableSampler();

    void calculate(Links *lattice);

    // Setters
    void setLatticeSize(int latticeSize);
    void setN(unsigned int *N);

    // Getters
    double getPlaquette() { return m_P; }
    double getTopologicalCharge() { return m_Q; }
    double getEnergyDensity() { return m_E; }
};

#endif // ObservableSampler_H
