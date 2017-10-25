#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"
#include "links.h"
#include "parallelization/indexorganiser.h"
#include <vector>

class Plaquette : public Correlator
{
private:
    int muIndex[4];
    int nuIndex[4];
    double multiplicationFactor;
    SU3 P;
public:
    Plaquette();
    ~Plaquette();
    double calculate(Links *lattice);
    void setLatticeSize(int latticeSize);

    inline void updateMuIndex(int mu) {
        for (int i = 0; i < 4; i++)
        {
            muIndex[i] = 0;
        }
        muIndex[mu] = 1;
    }
    inline void updateNuIndex(int nu) {
        for (int i = 0; i < 4; i++)
        {
            nuIndex[i] = 0;
        }
        nuIndex[nu] = 1;
    }
};

#endif // PLAQUETTE_H
