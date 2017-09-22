#ifndef WILSONGAUGEACTION_H
#define WILSONGAUGEACTION_H

#include "action.h"

class WilsonGaugeAction : public Action
{
private:
    // Lorentz indices arrays
    int *muIndex;
    int *nuIndex;
    // Action based constants
    double m_beta;
    double multiplicationFactor;
    SU3 staple;
    SU3 tr;
public:
    WilsonGaugeAction(double beta);
    ~WilsonGaugeAction();
    double getDeltaAction(Links *lattice, SU3 UPrime, int i, int j, int k, int l, int mu);
    void computeStaple(Links *lattice, int i, int j, int k, int l, int mu);
};

#endif // WILSONGAUGEACTION_H
