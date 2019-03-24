/**
 * \class WilsonGaugeAction
 *
 * \brief An implementation of the Wilson gauge action,
 *
 * \f{eqnarray*}{
 * S_G[U] = \frac{\beta}{3} \sum_{n\in\Lambda} \sum_{\mu<\nu} \mathrm{Re} \mathrm{tr} \big[ 1 - P_{\mu\nu}(n) \big].
 * \f}
 *
 * The different between this and LuscherAction is that in the
 * getActionDerivative we perform an explicit calculation of the derivative
 * in terms of its \f$\mathrm{SU}(3)\f$ generator.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
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
    SU3 m_staple, m_staple1, m_staple2;
    Lattice<SU3> m_latticeStaple,m_tempStaple1,m_tempStaple2;
    Lattice<double> m_tempDiag;
public:
    WilsonGaugeAction();
    ~WilsonGaugeAction();
    double getDeltaAction(SU3 U, SU3 UPrime);
    void computeStaple(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    Lattice<SU3> getActionDerivative(Lattice<SU3> *lattice, int mu);

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
