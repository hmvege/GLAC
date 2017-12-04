#ifndef WILSONGAUGEACTION_H
#define WILSONGAUGEACTION_H

#include "action.h"

class WilsonGaugeAction : public Action
{
private:
    // Lorentz indices arrays
    int m_muIndex[4];
    int m_nuIndex[4];
    // Action based constants
    double m_beta;
    double m_multiplicationFactor;
    SU3 m_staple, m_X, m_staple1, m_staple2;

    // FOR ALTERNATIVE DERIVATIVE METHOD!
    SU3 Q, C, Omega;
    double tempDiag;

public:
    WilsonGaugeAction();
    ~WilsonGaugeAction();
    double getDeltaAction(SU3 U, SU3 UPrime);
    void computeStaple(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    SU3 getActionDerivative(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);

    inline void updateMuIndex(int mu) {
        for (int i = 0; i < 4; i++)
        {
            m_muIndex[i] = 0;
        }
        m_muIndex[mu] = 1;
    }
    inline void updateNuIndex(int nu) {
        for (int i = 0; i < 4; i++)
        {
            m_nuIndex[i] = 0;
        }
        m_nuIndex[nu] = 1;
    }
};

#endif // WILSONGAUGEACTION_H
