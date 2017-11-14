#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"

class Plaquette : public Correlator
{
private:
    double m_multiplicationFactor;
    SU3 P;
    SU3 PTemp;
    const static std::string m_observableName;
public:
    Plaquette();
    ~Plaquette();
    void calculate(Links *lattice, int i);
    void calculate(SU3 *plaquetteStaples, int i);
    void setLatticeSize(int latticeSize);
    void setPlaquetteStaples(SU3 *plaquetteStaples);
};

#endif // PLAQUETTE_H
