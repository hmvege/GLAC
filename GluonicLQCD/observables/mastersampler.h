#ifndef MASTERSAMPLER_H
#define MASTERSAMPLER_H

#include "math/lattice.h"
#include "math/latticemath.h"

class MasterSampler
{
private:
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor;
    int m_latticeSize;
public:
    MasterSampler();
    ~MasterSampler() {}
    void calculate();
};

#endif // MASTERSAMPLER_H
