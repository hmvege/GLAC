#ifndef MASTERSAMPLER_H
#define MASTERSAMPLER_H

#include "math/lattice.h"
#include "math/latticemath.h"

class MasterSampler
{
private:
    double m_multiplicationFactor;
//    SU3 P;
//    SU3 PTemp;
public:
    MasterSampler();
    ~MasterSampler() {}
    void calculate();
};

#endif // MASTERSAMPLER_H
