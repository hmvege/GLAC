#ifndef WILSONGAUGEACTION_H
#define WILSONGAUGEACTION_H

#include "action.h"

class WilsonGaugeAction : public Action
{
private:
    // Lorentz indices arrays
//    int *muIndex;
//    int *nuIndex;
    int muIndex[4];
    int nuIndex[4];
    // Action based constants
    double m_beta;
    double multiplicationFactor;
    SU3 staple;
//    SU3 tr;
public:
    WilsonGaugeAction(double beta);
    ~WilsonGaugeAction();
    double getDeltaAction(Links *lattice, SU3 UPrime, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    void computeStaple(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);

    inline void updateMuIndex(int mu) {
        for (int i = 0; i < 4; i++)
        {
            muIndex[i] = 0;
        }
        muIndex[mu] = 1;
    }
    inline void updateNuIndex(int nu) {
        for (int i = 0; i < 4; i++)
        {
            nuIndex[i] = 0;
        }
        nuIndex[nu] = 1;
    }
};

#endif // WILSONGAUGEACTION_H
