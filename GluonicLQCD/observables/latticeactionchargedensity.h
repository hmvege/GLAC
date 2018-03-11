#ifndef LATTICEACTIONCHARGEDENSITY_H
#define LATTICEACTIONCHARGEDENSITY_H

#include "math/lattice.h"
#include "observables/observables.h"

class LatticeActionChargeDensity : public Correlator
{
private:
    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    double m_plaquette;
    Lattice <double> m_tempDiag, m_topCharge, m_energy;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
public:
    LatticeActionChargeDensity(bool flow);
    ~LatticeActionChargeDensity();
    void calculate(Lattice<SU3> * lattice, int iObs);
    void initializeObservableStorer(bool storeFlowObservable);

    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(int iFlow);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(int iObs);
    void printStatistics();
    std::vector<double> getObservablesVector(int iObs);
    void copyObservable(int iObs, std::vector<double> obs);
};

#endif // LATTICEACTIONCHARGEDENSITY_H
