#ifndef OBSERVABLESAMPLER_H
#define OBSERVABLESAMPLER_H

#include "correlator.h"
#include "plaquette.h"
#include "clover.h"
#include "energydensity.h"
#include "topologicalcharge.h"

class ObservableSampler : public Correlator
{
private:
    unsigned int m_N[4];
    int m_latticeSize;
    double m_a;
    double m_P, m_E, m_Q;
    std::vector<int> m_position;

    // All observables in program
    Clover *m_clover = nullptr;
    Plaquette *m_plaquette = nullptr;
    EnergyDensity *m_energyDensity = nullptr;
    TopologicalCharge *m_topologicalCharge = nullptr;
public:
    ObservableSampler(bool storeFlowObservable);
    ~ObservableSampler();

    void calculate(Links *lattice, int iObs);

    void writeStatisticsToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(int iConfig);
    void runStatistics();

    // Printers
    void printObservable(int iObs);
    void printHeader();

    // Setters
    void reset();
    void setLatticeSize(int latticeSize);
    void setN(unsigned int *N);

    // Getters
    double getObservable(int iObs);
    double getPlaquette() { return m_P; }
    double getTopologicalCharge() { return m_Q; }
    double getEnergyDensity() { return m_E; }
};

#endif // ObservableSampler_H
