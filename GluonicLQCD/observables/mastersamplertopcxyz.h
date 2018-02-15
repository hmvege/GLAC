#ifndef MASTERSAMPLERTOPCXYZ_H
#define MASTERSAMPLERTOPCXYZ_H

#include "math/lattice.h"
#include "observables/observables.h"

class MasterSamplerTopcXYZ : public Correlator
{
private:
    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    double m_topCharge, m_energy, m_plaquette;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Container for the topc xyz observable
//    std::vector<std::vector<double>> m_topcxyz;
    std::vector<double> m_tempTopcT;
    double * m_tempTopctArray;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
    ObservableStorer * m_topctObservable = nullptr;
public:
    MasterSamplerTopcXYZ(bool flow);
    ~MasterSamplerTopcXYZ();
    void calculate(Lattice<SU3> * lattice, int iObs);
    void storeFlow(bool storeFlowObservable);

    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(int configNumber);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(int iObs);
    void printStatistics();
    std::vector<double> getObservablesVector(int iObs);
    void copyObservable(int iObs, std::vector<double> obs);
};

#endif // MASTERSAMPLERTOPCXYZ_H
