#include "flow.h"
#include "links.h"
#include "parallelization/indexorganiser.h"
#include "functions.h"

// FOR TESTING IF FLOW PRESERVES SU3
#include "testsuite.h"

/*
 * A class for performing Wilson flow on a gauge field configuration.
 */

Flow::Flow(unsigned int *N, double beta)
{
    for (int i = 0; i < 4; i++) m_N[i] = N[i];
    m_subLatticeSize = 1;
    for (int i = 0; i < 4; i++) m_subLatticeSize *= m_N[i];
    f0.identity();
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
     * Performs a NFlows of flow on the lattice.
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
//                        W[0] = lattice[m_Index->getIndex(x,y,z,t)].U[mu]; // V should be the previous flowed point!
                        m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_S->getActionDerivative(lattice,x,y,z,t,mu)*m_epsilon);

//                        m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_S->getActionDerivative(lattice,x,y,z,t,mu));
//                        SU3 U = exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]);
//                        SU3 UInv = exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*(-1.0));
//                        SU3 I = U*U.inv();
//                        I.print();exit(1);
//                        exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]).print();exit(1);
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
                        std::cout<<std::endl;
                        (m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]).printMachine();
                        std::cout<<std::endl;
                        exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]).printMachine();
                        std::cout<<std::endl;
                        exponentiate2(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]).printMachine();
                        std::cout<<std::endl;
                        exit(1);
                        lattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]*0.25)*lattice[m_Index->getIndex(x,y,z,t)].U[mu]);

                        (exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu])*exponentiate(m_tempLattice[m_Index->getIndex(x,y,z,t)].U[mu]).inv()).print();exit(1);
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
     * Takes the exponential of an Hermitian matrix.
     */
    f0.identity();

//    // Makes Q hermitian. MOVE TO FUNCTION?
//    for (int i = 0; i < 9; i++) {
//        temp = Q.mat[2*i];
//        Q.mat[2*i] = Q.mat[2*i+1];
//        Q.mat[2*i+1] = -temp;
//    }

    QSquared = Q*Q;
    QCubed = Q*QSquared;


//    std::cout << "Q: " << std::endl;
//    Q.print();
//    std::cout << "QSquared: " << std::endl;
//    QSquared.print();
//    std::cout << "QCubed: " << std::endl;
//    QCubed.print();

    c0 = 0.3333333333333333*(QCubed.mat[0] + QCubed.mat[8] + QCubed.mat[16]);
    c1 = 0.5*(QSquared.mat[0] + QSquared.mat[8] + QSquared.mat[16]);

//    std::cout << "c0 = " << c0 << std::endl;
//    std::cout << "c1 = " << c1 << std::endl;

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
        xi0 = 1 - 0.16666666666666666*ww*(1 - 0.05*ww*(1 - 0.023809523809523808*ww)); // 1/6.0 and 1/42.0
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
        for (int i = 0; i < 3; i++) {
            f[i] = f[i].c()*pow(-1,double(i));
        }
    }

    // Testing if we uphold the identity Q^3 - c1*Q - c0*I = 0 =============
//    SU3 X = QCubed - Q*c1 - I*c0;
//    for (int i = 0; i < 18; i++) {
//        if (X.mat[i] > 1e-16) {
//            std::cout << "ERROR: \nX: " << std::endl;
//            X.print();
//            std::cout << "EXITS" << std::endl;
//            exit(1);
//        }
//    }
    // =====================================================================

    f0.setComplex(f[0],0);
    f0.setComplex(f[0],8);
    f0.setComplex(f[0],16);
    return f0 + Q*f[1] + QSquared*f[1];
}

SU3 Flow::exponentiate2(SU3 Q)
{
    Q.makeHermitian(); // Makes anti-hermitian
    SU3 U1,U2,U3,Y1,Y2,Y3,I;
    I.identity();
    Y1.zeros();
    Y2.zeros();
    Y3.zeros();
    complex x1,x2,x3;
    x1 = complex(Q.mat[0] - Q.mat[8],Q.mat[1] - Q.mat[9]) * 0.3333333333333333;
    Y1.setComplex(x1,0);
    Y1.setComplex(complex(-x1.re(),-x1.im()),4);
    Y1.mat[2] = Q.mat[2];
    Y1.mat[3] = Q.mat[3];
    Y1.mat[6] = Q.mat[6];
    Y1.mat[7] = Q.mat[7];

    x2 = complex(Q.mat[0] - Q.mat[16],Q.mat[1] - Q.mat[17]) * 0.3333333333333333;
    Y2.setComplex(x2,0);
    Y2.setComplex(complex(-x2.re(),-x2.im()),8);
    Y2.mat[4] = Q.mat[4];
    Y2.mat[5] = Q.mat[5];
    Y2.mat[12] = Q.mat[12];
    Y2.mat[13] = Q.mat[13];

    x3 = complex(Q.mat[8] - Q.mat[16],Q.mat[9] - Q.mat[17]) * 0.3333333333333333;
    Y3.setComplex(x3,4);
    Y3.setComplex(complex(-x3.re(),-x3.im()),8);
    Y3.mat[10] = Q.mat[10];
    Y3.mat[11] = Q.mat[11];
    Y3.mat[14] = Q.mat[14];
    Y3.mat[15] = Q.mat[15];

//    * 0 1 2    0  1   2  3    4  5
//    * 3 4 5 =  6  7   8  9   10 11
//    * 6 7 8   12 13  14 15   16 17

//    U1 = (I + Y1*0.25) * inverseY1;
//    U1 = (I + Y2*0.25) * inverseY2;
//    U1 = (I + Y3*0.5) * inverseY3;

    U1 *= U2;
    U1 *= U3;
    U1 *= U2;
    U1 *= U1;
    return U1;
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
