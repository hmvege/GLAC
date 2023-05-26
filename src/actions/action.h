/*!
 * \class Action
 *
 * \brief The Action class
 *
 * The base action class for all other actions to be derived upon.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef ACTION_H
#define ACTION_H

#include "math/lattice.h"
#include <vector>

using std::cout;
using std::endl;

class Action
{
protected:
    /*!
     * \brief m_N Lattice dimension array
     */
    std::vector<unsigned int> m_N;
    //! Array for holding the position in lattice.
    /*!
     * \brief m_position is used for handling the shift-method in parallelization.
     */
    std::vector<int> m_position;
public:

    //! \brief Action base constructor
    Action();

    /*!
     * \brief ~Action destructor. Nothing to de-allocate in base class
     */
    virtual ~Action();

    /*!
     * \brief getDeltaAction computes the change in action. Staple must have already have been calculated.
     * \param U Old link
     * \param UPrime Updated link
     */
    virtual double getDeltaAction(SU3 &U, SU3 &UPrime) = 0;

    /*!
     * \brief computeStaple computes the staple at given poisition.
     * \param lattice a pointer of four lattice objects, one for each lorentz index.
     * \param i spatial \f$x\f$ position.
     * \param j spatial \f$y\f$ position.
     * \param k spatial \f$z\f$ position.
     * \param l temporal \f$t\f$ position.
     * \param mu lorentz index, \f$\mu\f$.
     * \return Returns change in action \f$\Delta S\f$
     */
    virtual void computeStaple(Lattice<SU3> *lattice, const int i, const int j, const int k, const int l, const int mu) = 0;

    /*!
     * \brief getActionDerivative, computes the derivative of the lattice in given direction.
     * \param lattice a pointer of four lattice objects, one for each lorentz index.
     * \param mu lorentz index, \f$\mu\f$.
     */
    virtual Lattice<SU3> getActionDerivative(Lattice<SU3> * lattice, const int mu) = 0;
};



#endif // ACTION_H
