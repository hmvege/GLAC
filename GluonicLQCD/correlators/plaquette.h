#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"
#include "links.h"

class Plaquette : public Correlator
{
private:
    int *muIndex;
    int *nuIndex;
public:
    Plaquette(int N, int N_T);
    ~Plaquette();
    double calculate(Links *lattice);
};

#endif // PLAQUETTE_H
