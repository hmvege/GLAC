#ifndef WILSONGAUGEACTION_H
#define WILSONGAUGEACTION_H

#include "action.h"

class WilsonGaugeAction : public Action
{
private:
    double m_beta;
public:
    WilsonGaugeAction(int latticeSize, double beta);
    ~WilsonGaugeAction();
    double getDeltaAction(Links *lattice, SU3 U, int i, int j, int k, int l, int mu);
};

#endif // WILSONGAUGEACTION_H
