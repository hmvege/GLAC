#ifndef MASTERSAMPLER_H
#define MASTERSAMPLER_H

#include "math/lattice.h"
#include "observables/observables.h"

class MasterSampler : public Correlator
{
private:
    int mu;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    double m_topCharge, m_energy, m_plaquette;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
public:
    MasterSampler(bool flow);
    ~MasterSampler() {}
    void calculate(Lattice<SU3> * lattice, int iObs);
    void storeFlow(bool storeFlowObservable);

    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(int iFlow);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(int iObs);
    void printStatistics();
};

#endif // MASTERSAMPLER_H
