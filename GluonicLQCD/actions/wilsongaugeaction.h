#ifndef WILSONGAUGEACTION_H
#define WILSONGAUGEACTION_H

#include "action.h"

class WilsonGaugeAction : public Action
{
private:
    double m_beta;
    SU3 m_staple;
    double multiplicationFactor;
public:
    WilsonGaugeAction(int N, int N_T, double beta);
    ~WilsonGaugeAction();
    double getDeltaAction(Links *lattice, SU3 UPrime, int i, int j, int k, int l, int mu);
    void computeStaple(Links *lattice, int i, int j, int k, int l, int mu);
};

#endif // WILSONGAUGEACTION_H
