#ifndef TOPCHARGEANDENERGY_H
#define TOPCHARGEANDENERGY_H

#include "energydensity.h"
#include "topchargeandenergy.h"
#include "correlator.h"

class TopChargeAndEnergy : public Correlator
{
public:
    TopChargeAndEnergy();
    ~TopChargeAndEnergy();

    double calculate(Links *lattice);
};

#endif // TOPCHARGEANDENERGY_H
