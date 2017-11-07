#ifndef CLOVER_H
#define CLOVER_H

#include "correlator.h"

class Clover : public Correlator
{
private:
    SU3 U1,U2,U3,U4,U2Temp,U3Temp;
    int m_cloverOverCounter = 0;
    int m_plaquetteOverCounter = 0;
public:
    Clover();
    ~Clover();

    SU3 m_clovers[12];
    SU3 m_plaquettes[6];
    void calculateClover(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l);
};

#endif // CLOVER_H
