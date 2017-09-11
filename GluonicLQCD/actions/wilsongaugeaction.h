#ifndef WILSONGAUGEACTION_H
#define WILSONGAUGEACTION_H

#include "action.h"
#include <vector>

class WilsonGaugeAction : public Action
{
private:
    double m_beta;
    double multiplicationFactor;
    SU3 m_staple;
    std::vector<int> indexes;
public:
    WilsonGaugeAction(double beta);
    ~WilsonGaugeAction();
    double getDeltaAction(Links *lattice, SU3 UPrime, int i, int j, int k, int l, int mu);
    void computeStaple(Links *lattice, int i, int j, int k, int l, int mu);
};

#endif // WILSONGAUGEACTION_H
