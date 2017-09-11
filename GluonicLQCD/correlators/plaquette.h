#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"
#include "links.h"
#include "parallelization/indexorganiser.h"
#include <vector>

class Plaquette : public Correlator
{
private:
    int *muIndex;
    int *nuIndex;
    double multiplicationFactor;
public:
//    Plaquette(int N, int N_T);
    Plaquette();
    ~Plaquette();
    double calculate(Links *lattice);
    void setLatticeSize(int latticeSize);
};

#endif // PLAQUETTE_H
