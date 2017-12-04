#ifndef MASTERSAMPLER_H
#define MASTERSAMPLER_H

#include "math/lattice.h"
#include "observables/observables.h"

class MasterSampler : public Correlator
{
private:
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
//    Lattice<SU3> clov1(dim), clov2(dim), U2Temp(dim), U3Temp(dim), Temp1(dim);
    double m_topCharge, m_energy, m_plaquette;
    Lattice <double> m_tempDiag;
    int mu;
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

public:
    MasterSampler(bool flow);
    ~MasterSampler() {}
    void calculate(Lattice<SU3> * lattice, int iObs);
};

#endif // MASTERSAMPLER_H
