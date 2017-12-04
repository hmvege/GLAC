#ifndef MASTERSAMPLER_H
#define MASTERSAMPLER_H

#include "math/lattice.h"


class MasterSampler
{
private:
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    int m_latticeSize;

//    Lattice<SU3> clov1(dim), clov2(dim), U2Temp(dim), U3Temp(dim), Temp1(dim);
public:
    MasterSampler();
    ~MasterSampler() {}
    void calculate();
};

#endif // MASTERSAMPLER_H
