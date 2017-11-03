#include "flow.h"
#include "links.h"
#include "parallelization/indexorganiser.h"
#include "functions.h"
#include <mpi.h>

// FOR TESTING IF FLOW PRESERVES SU3
#include "testsuite.h"
#include <iostream>

using std::cout;
using std::endl;

/*
 * A class for performing Wilson flow on a gauge field configuration.
 */

Flow::Flow(unsigned int *N, double beta, int numprocs, int processRank)
{
    m_numprocs = numprocs;
    m_processRank = processRank;
    for (int i = 0; i < 4; i++) m_N[i] = N[i];
    m_subLatticeSize = 1;
    for (int i = 0; i < 4; i++) m_subLatticeSize *= m_N[i];
    f0.identity();
    I.identity();
    m_beta = beta;
    m_tempLattice = new Links[m_subLatticeSize];
}

Flow::~Flow()
{
    delete [] m_tempLattice;
}

void Flow::flowGaugeField(int NFlows, Links *lattice)
{
    /*
     * Performs a NFlows of flow on the lattice. REMOVE THIS SINCE IT IS REDUNDANT! ALWAYS SAMPLING THE FLOW ITERATION!
     */
    for (int i = 0; i < NFlows; i++) {
        runFlow(lattice);
    }

//    // Tests if flow preserves the SU3 matrix properties ===================
//    TestSuite test;
//    for (unsigned int i = 0; i < m_subLatticeSize; i++) {
//        for (unsigned int mu = 0; mu < 4; mu++) {
//            test.testMatrix(lattice[i].U[mu],true);
//            exit(1);
//        }
//    }
//    // =====================================================================

}

void Flow::runFlow(Links *lattice)
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
                        m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_S->getActionDerivative(lattice,x,y,z,t,mu)*m_epsilon);
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
//                        if (m_processRank == 0) {
//                            (m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]).printMachine();
//                            cout<<"MORNINGSTAR METHOD: "<<endl;
//                            exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*0.25).printMachine();
//                            cout<<"LUSCHER METHOD: "<<endl;
//                            exponentiate2(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*0.25).printMachine();
//                            cout<<"TAYLOR EXPANSION: "<<endl;
//                            exponentiate3(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*0.25).printMachine();
//                        }
//                        MPI_Finalize();exit(1);
                        lattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*0.25)*lattice[m_Index->getIndex(x,y,z,t)].U[mu]);
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
                        m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_S->getActionDerivative(lattice,x,y,z,t,mu)*m_epsilon*0.8888888888888888 - m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*0.4722222222222222);
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
                        lattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu])*lattice[m_Index->getIndex(x,y,z,t)].U[mu]);
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
                        m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_S->getActionDerivative(lattice,x,y,z,t,mu)*m_epsilon*0.75 - m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]);
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
                        lattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu])*lattice[m_Index->getIndex(x,y,z,t)].U[mu]);
                    }
                }
            }
        }
    }
//    lattice[m_Index->getIndex(0,0,0,0)].U[0].print();
}

SU3 Flow::exponentiate(SU3 Q)
{
    /*
     * Takes the exponential of an Hermitian matrix using the Morningstar method.
     */
    f0.zeros();
    Q.makeHermitian();

    QSquared = Q*Q;
    QCubed = Q*QSquared;

    c0 = 0.3333333333333333*(QCubed.mat[0] + QCubed.mat[8] + QCubed.mat[16]);
    c1 = 0.5*(QSquared.mat[0] + QSquared.mat[8] + QSquared.mat[16]);

    c0max = 0.6666666666666666 * c1 * sqrt(c1 * 0.3333333333333333);
    theta = acos(fabs(c0)/c0max);
    u = sqrt(0.3333333333333333*c1) * cos(0.3333333333333333 * theta);
    w = sqrt(c1) * sin(0.3333333333333333 * theta);

    // Sets shortenings
    uu = u*u;
    ww = w*w;
    cosu = cos(u);
    sinu = sin(u);
    cosw = cos(w);
    sinw = sin(w);
    cos2u = cos(2*u);
    sin2u = sin(2*u);

    // Find xi(w)
    if (fabs(w) > 0.05) {
        xi0 = sinw/w;
    } else {
        xi0 = 1 - ww*(1 - 0.05*ww*(1 - ww/42.0))/6.0; // 1/6.0 and 1/42.0
    }

    // Sets h
    h[0].setRe(8*uu*cosu*cosw + 2*u*xi0*(3*uu + ww)*sinu + (uu - ww)*cos2u);
    h[0].setIm(-8*uu*sinu*cosw + 2*u*xi0*(3*uu + ww)*cosu + (uu - ww)*sin2u);
    h[1].setRe(-2*u*cosu*cosw + 2*u*cos2u + xi0*(3*uu - ww)*sinu);
    h[1].setIm(2*u*sinu*cosw + 2*u*sin2u + xi0*(3*uu - ww)*cosu);
    h[2].setRe(-3*u*xi0*sinu - cosu*cosw + cos2u);
    h[2].setIm(-3*u*xi0*cosu + sinu*cosw + sin2u);

    // Sets f
    for (int i = 0; i < 3; i++) {
        f[i] = h[i] / (9*uu - ww);
    }

    // Checks for negative c0 coefficient
    if (c0 < 0) {
        f[0].conjugate();
        f[1] = -f[1].c();
        f[2].conjugate();
    }

    // Sets the first matrix, I*f0
    f0.setComplex(f[0],0);
    f0.setComplex(f[0],8);
    f0.setComplex(f[0],16);
    return f0 + Q*f[1] + QSquared*f[2];
}

SU3 Flow::exponentiate2(SU3 Q)
{
    /*
     * Exponentiation using the Luscher method.
     */
    // Ensures the U's are at zero for later filling
    U1.zeros();
    U2.zeros();
    U3.zeros();

    // Sets elements of the Y matrices
    x1 = (Q.get(0,0) - Q.get(1,1))*0.3333333333333333;
    x2 = (Q.get(0,0) - Q.get(2,2))*0.3333333333333333;
    x3 = (Q.get(1,1) - Q.get(2,2))*0.3333333333333333;

    // Sets often used factors
    X1221X = Q.get(0,1)*Q.get(1,0);
    X1331X = Q.get(0,2)*Q.get(2,0);
    X2332X = Q.get(1,2)*Q.get(2,1);

    // Sets complex denominators
    div1 = (X1221X + x1*x1)*0.0625 - 1.0;
    div2 = (X1331X + x2*x2)*0.0625 - 1.0;
    div3 = (X2332X + x3*x3)*0.25 - 1.0;

    // Setting U1 matrix
    sqrdFactor = x1*0.25 - 1.0;
    U1.setComplex((X1221X*0.0625 + x1*x1*0.0625 + x1*0.5 + 1.0)/div1*(-1.0),0);
    U1.setComplex((Q.get(0,1)*(-0.5))/div1,2);
    U1.setComplex((Q.get(1,0)*(-0.5))/div1,6);
    U1.setComplex((X1221X*0.0625 + sqrdFactor*sqrdFactor)/div1*(-1.0),8);
    U1.mat[16] = 1.0;

    // Setting U2 matrix
    sqrdFactor = x2*0.25 - 1.0;
    U2.setComplex((X1331X*0.0625 + x2*x2*0.0625 + x2*0.5 + 1.0)/div2*(-1),0);
    U2.setComplex((Q.get(0,2)*(-0.5))/div2,4);
    U2.mat[8] = 1.0;
    U2.setComplex((Q.get(2,0)*(-0.5))/div2,12);
    U2.setComplex((X1331X*0.0625 + sqrdFactor*sqrdFactor)/div2*(-1.0),16);

    // Setting U3 matrix
    sqrdFactor = x3*0.5 - 1.0;
    U3.mat[0] = 1.0;
    U3.setComplex((X2332X*0.25 + x3*x3*0.25 + x3 + 1.0)/div3*(-1.0),8);
    U3.setComplex(Q.get(1,2)/div3*(-1.0),10);
    U3.setComplex(Q.get(2,1)/div3*(-1.0),14);
    U3.setComplex((X2332X*0.25 + sqrdFactor*sqrdFactor)/div3*(-1.0),16);

    E = U1;
    E *= U2;
    E *= U3;
    E *= U2;
    E *= U1;
    return E;
}

SU3 Flow::exponentiate3(SU3 Q)
{
    /*
     * Exponentiate using regular Taylor expansion.
     */
    QSquared = Q;
    QSquared *= Q;
    QCubed = QSquared;
    QCubed *= Q;
    QQuartic = QCubed;
    QQuartic *= Q;
    return I + Q + QSquared*0.5 + QCubed/6.0 + QQuartic/24.0;
}

void Flow::setIndexHandler(IndexOrganiser *Index)
{
    /*
     * Sets the index handler that we are using. Must be set before performing the flow, or nonsense will ensue.
     */
    m_Index = Index;
}

void Flow::setAction(Action *S)
{
    /*
     * Sets the action class of the Flow.
     */
    m_S = S;
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
                        lattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]);
                    }
                }
            }
        }
    }
}
