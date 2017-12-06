#include "flow.h"
#include "parallelization/index.h"
#include "config/parameters.h"

/*
 * A class for performing Wilson flow on a gauge field configuration.
 */

Flow::Flow(Action *S)
{
    m_N = Parameters::getN();
    m_epsilon = Parameters::getFlowEpsilon();
    m_subLatticeSize = Parameters::getSubLatticeSize();
    m_tempLattice = new Lattice<SU3>[4];
    m_tempExpLattice.allocate(m_N);
    for (int mu = 0; mu < 4; mu++) {
        m_tempLattice[mu].allocate(m_N);
    }
    m_S = S;
    setSU3ExpFunc();
}

Flow::~Flow()
{
    delete [] m_tempLattice;
}

void Flow::flowField(Lattice<SU3> *lattice)
{
    /*
     * Performs a single flow on the lattice.
     */
    // W0 is simply just the original lattice times epsilon
    // Sets Z0 in temporary lattice
    for (unsigned int mu = 0; mu < 4; mu++) {
        m_tempLattice[mu] = m_S->getActionDerivative(lattice,mu);
    }
    // Sets W1 in main lattice
    for (unsigned int mu = 0; mu < 4; mu++) {
        lattice[mu] = matrixExp(m_tempLattice[mu]*(m_epsilon*0.25))*lattice[mu];
    }
    // Sets "Z1" in temporary lattice
    for (unsigned int mu = 0; mu < 4; mu++) {
        m_tempLattice[mu] = m_S->getActionDerivative(lattice,mu)*(m_epsilon*0.8888888888888888) - m_tempLattice[mu]*(m_epsilon*0.4722222222222222);
    }
    // Sets W2 in main lattice
    for (unsigned int mu = 0; mu < 4; mu++) {
        lattice[mu] = matrixExp(m_tempLattice[mu])*lattice[mu];
    }
    // Sets "Z2" in temporary lattice
    for (unsigned int mu = 0; mu < 4; mu++) {
        m_tempLattice[mu] = m_S->getActionDerivative(lattice,mu)*(0.75*m_epsilon) - m_tempLattice[mu];
    }
    // Sets V_{t+1} in main lattice
    for (unsigned int mu = 0; mu < 4; mu++) {
        lattice[mu] = matrixExp(m_tempLattice[mu])*lattice[mu];
    }
}

void Flow::setSU3ExpFunc()
{
    /*
     * Sets the prefered method of exponentiation.
     * Does nothing if name is not recognized.
     * If Morningstar, do nothing as that is default.
     */
    if (Parameters::getExpFuncName() == "taylor2") { // ENSURE WE DONT GO OUT OF SCOPE!!
        m_SU3ExpFunc = new Taylor2Exp;
    } else if (Parameters::getExpFuncName() == "taylor4") {
        m_SU3ExpFunc = new Taylor4Exp;
    } else if (Parameters::getExpFuncName() == "luscher") {
        m_SU3ExpFunc = new ExpLuscher;
    } else if (Parameters::getExpFuncName() == "morningstar") {
        m_SU3ExpFunc = new SU3Exp;
    } else {
        printf("SU3 exp. func. %s not recognized",Parameters::getExpFuncName().c_str());
    }
}

inline Lattice<SU3> Flow::matrixExp(const Lattice<SU3> &lattice)
{
    for (unsigned int iSite = 0; iSite < lattice.m_latticeSize; iSite++) {
        m_tempExpLattice.m_sites[iSite] = m_SU3ExpFunc->exp(lattice.m_sites[iSite]);
    }
    return m_tempExpLattice;
}

inline Lattice<SU3> Flow::matrixExp(Lattice<SU3> &&lattice)
{
    for (unsigned int iSite = 0; iSite < lattice.m_latticeSize; iSite++) {
        m_tempExpLattice.m_sites[iSite] = m_SU3ExpFunc->exp(lattice.m_sites[iSite]);
    }
    return m_tempExpLattice;
}
