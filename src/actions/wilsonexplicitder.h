/**
 * \class WilsonExplicitDer
 *
 * \brief An implementation of the Wilson gauge action,
 *
 * \f{eqnarray*}{
 * S_G[U] = \frac{\beta}{3} \sum_{n\in\Lambda} \sum_{\mu<\nu} \mathrm{Re} \mathrm{tr} \big[ 1 - P_{\mu\nu}(n) \big].
 * \f}
 *
 * The different between this and WilsonGaugeAction is that in the
 * getActionDerivative we perform an explicit calculation of the derivative
 * in terms of its \f$\mathrm{SU}(3)\f$ generator.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef WILSONEXPLICITDER_H
#define WILSONEXPLICITDER_H

#include "action.h"

class WilsonExplicitDer : public Action
{
private:
    /*!
     * \brief m_muIndex Lorentz indices array.
     *
     * Needed for loop over the two Lorentz indices in the field strength tensor
     */
    int m_muIndex[4];

    /*!
     * \brief m_nuIndex Lorentz indices array.
     *
     * Needed for loop over the two Lorentz indices in the field strength tensor
     */
    int m_nuIndex[4];

    /*!
     * \brief m_beta the strong coupling constant.
     *
     * \f$\beta\f$ determines the strength of interaction and the lattice spacing(and more).
     */
    double m_beta;

    /*!
     * \brief m_multiplicationFactor to be factored in at the end.
     *
     * m_multiplicationFactor is setup during construction. Saves a few flops.
     */
    double m_multiplicationFactor;

    /*!
     * \brief m_staple, m_staple1, m_staple2 is needed for computing the staple.
     */
    SU3 m_staple, m_staple1, m_staple2;

    /*!
     * \brief m_latticeStaple
     */
    Lattice<SU3> m_latticeStaple,m_tempStaple1,m_tempStaple2;

public:
    WilsonExplicitDer();
    ~WilsonExplicitDer();
    double getDeltaAction(SU3 &U, SU3 &UPrime);
    void computeStaple(Lattice<SU3> *lattice, const int i, const int j, const int k, const int l, const int mu);
    Lattice<SU3> getActionDerivative(Lattice<SU3> *lattice, const int mu);

    /*!
     * \brief updateMuIndex updates the m_muIndex based on mu.
     * \param mu direction which will be set to 1. All others set to zero.
     */
    inline void updateMuIndex(int mu) {
        for (int i = 0; i < 4; i++)
        {
            m_muIndex[i] = 0;
        }
        m_muIndex[mu] = 1;
    }

    /*!
     * \brief updateNuIndex updates the m_nuIndex based on nu.
     * \param nu direction which will be set to 1. All others set to zero.
     */
    inline void updateNuIndex(int nu) {
        for (int i = 0; i < 4; i++)
        {
            m_nuIndex[i] = 0;
        }
        m_nuIndex[nu] = 1;
    }
};

#endif // WILSONEXPLICITDER_H
