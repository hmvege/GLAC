#include "flow.h"

/*
 * A class for performing Wilson flow on a gauge field configuration.
 */

Flow::Flow()
{
    Parameters::getN(m_N);
    m_epsilon = Parameters::getFlowEpsilon();
    m_subLatticeSize = Parameters::getSubLatticeSize();
    m_tempLattice = new Links[m_subLatticeSize];
}

Flow::Flow(Action *S)
{
    Parameters::getN(m_N);
    m_epsilon = Parameters::getFlowEpsilon();
    m_subLatticeSize = Parameters::getSubLatticeSize();
    m_tempLattice = new Links[m_subLatticeSize];
    m_S = S;
}

Flow::~Flow()
{
    delete [] m_tempLattice;
}

void Flow::flowField(Links *lattice)
{
    /*
     * Performs a single flow on the lattice.
     */
    // W0 is simply just the original lattice times epsilon
    // Sets Z0 in temporary lattice
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {

                        m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_S->getActionDerivative(lattice,x,y,z,t,mu);
                    }
                }
            }
        }
    }
    // Sets W1 in main lattice
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_SU3ExpFunc->exp(m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu]*m_epsilon*0.25)*lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu];
                    }
                }
            }
        }
    }
    // Sets "Z1" in temporary lattice
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_S->getActionDerivative(lattice,x,y,z,t,mu)*m_epsilon*0.8888888888888888 - m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu]*m_epsilon*0.4722222222222222;
                    }
                }
            }
        }
    }
    // Sets W2 in main lattice
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_SU3ExpFunc->exp(m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu])*lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu];
                    }
                }
            }
        }
    }
    // Sets "Z2" in temporary lattice
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_S->getActionDerivative(lattice,x,y,z,t,mu)*0.75*m_epsilon - m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu];
                    }
                }
            }
        }
    }
    // Sets V_{t+1} in main lattice
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_SU3ExpFunc->exp(m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu])*lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu];
                    }
                }
            }
        }
    }
}

void Flow::setAction(Action *S)
{
    /*
     * Sets the action class of the Flow.
     */
    m_S = S;
}

void Flow::setSU3ExpFunc(SU3Exp *SU3ExpFunc)
{
    /*
     * Sets the prefered method of exponentiation.
     */
    m_SU3ExpFunc = SU3ExpFunc;
}

inline void Flow::updateLattice(Links *lattice)
{
    /*
     * Updates the lattice new values from m_tempLattice.
     */
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_tempLattice[Parallel::Index::getIndex(x,y,z,t)].U[mu];
                    }
                }
            }
        }
    }
}
