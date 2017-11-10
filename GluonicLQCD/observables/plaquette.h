#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"

class Plaquette : public Correlator
{
private:
    double m_multiplicationFactor;
    SU3 P;
    SU3 PTemp;
public:
    Plaquette();
    ~Plaquette();
    double calculate(Links *lattice);
    double calculate(SU3 *plaquetteStaples);
    void setLatticeSize(int latticeSize);
    void setPlaquetteStaples(SU3 *plaquetteStaples);
};

#endif // PLAQUETTE_H
