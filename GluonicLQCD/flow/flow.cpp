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
    I.identity();
    m_beta = beta;
    m_updatedLattice = new Links[m_subLatticeSize];
}

Flow::~Flow()
{
    delete [] m_updatedLattice;
}

void Flow::flowGaugeField(int NFlows, Links *lattice)
{
    /*
     * Performs a NFlows of flow on the lattice.
     */
    for (int i = 0; i < NFlows; i++) {
        runFlow(lattice);
    }

    // Tests if flow preserves the SU3 matrix properties
    TestSuite test;
    for (unsigned int i = 0; i < m_subLatticeSize; i++) {
        for (unsigned int mu = 0; mu < 4; mu++) {
            test.testMatrix(lattice[i].U[mu],true);
            exit(1);
        }

    }
}

void Flow::runFlow(Links *lattice)
{
    /*
     * Performs a single flow on the lattice.
     */
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        smearLink(lattice,x,y,z,t,mu);
                    }
                }
            }
        }
    }
    updateLattice(lattice);
}

SU3 Flow::exponentiate(SU3 Q)
{
    /*
     * Takes the exponential of an Hermitian matrix.
     */
    // Makes Q hermitian. MOVE TO FUNCTION?
    I.identity();

    double temp = 0;
    for (int i = 0; i < 9; i++) {
        temp = Q.mat[2*i];
        Q.mat[2*i] = Q.mat[2*i+1];
        Q.mat[2*i+1] = -temp;
    }

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
//    h[0].setRe((u*u - w*w)*cos(2*w) + 8*u*u*cos(u)*cos(w) + 2*u*(3*u*u + w*w)*sin(w)*xi0); // SIMPLIFY w*w, sin(w), cos(w), cos(2*u), sin(2*u), 3*u*u - w*w
//    h[0].setIm((u*u - w*w)*sin(2*u) + cos(u)*2*u*(3*u*u + w*w)*xi0 - sin(w)*8*u*u*cos(w));
//    h[1].setRe(2*u*cos(2*u) - cos(u)*2*u*cos(w) - sin(u)*(3*u*u - w*w)*xi0);
//    h[1].setIm(2*u*sin(2*u) + cos(u)*(3*u*u - w*w)*xi0 - sin(u)*2*u*cos(w));
//    h[2].setRe(cos(2*u) - cos(u)*cos(w) - sin(u)*3*u*xi0);
//    h[2].setIm(sin(2*u) - cos(u)*3*u*xi0 + sin(u)*cos(w));

    // Sets f
    for (int i = 0; i < 3; i++) {
        f[i] = h[i] / (9*uu - ww);
        if (c0 < 0) {
            f[0] = f[0].c();
            f[1] = f[1].c()*(-1);
            f[2] = f[2].c();
        }
    }

    // Testing if we uphold the identity Q^3 - c1*Q - c0*I = 0
    SU3 X = QCubed - Q*c1 - I*c0;
    for (int i = 0; i < 18; i++) {
        if (X.mat[i] > 1e-16) {
            std::cout << "ERROR: \nX: " << std::endl;
            X.print();
            std::cout << "EXITS" << std::endl;
            exit(1);
        }
    }
//    exit(1);
//    I.setComplex(f[0],0);
//    I.setComplex(f[0],8);
//    I.setComplex(f[0],16);
    return I*f[0] + Q*f[1] + QSquared*f[1];
}

void Flow::smearLink(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    /*
     * Smears a link according to what that has been outlined in Luschers paper on Wilson flow.
     */
    // Sets first RK3 constant, W0
    W[0] = lattice[m_Index->getIndex(i,j,k,l)].U[mu]; // V should be the previous flowed point!
    // Finds Z0
    Z[0] = m_S->getActionDerivative(lattice,W[0],i,j,k,l,mu) * m_epsilon;
    // Sets second RK3 constant, W1
    W[1] = exponentiate(Z[0] * 0.25) * W[0];
    // Finds Z1
    Z[1] = m_S->getActionDerivative(lattice,W[1],i,j,k,l,mu) * m_epsilon;
    // Sets third RK3 constant, W2
    W[2] = exponentiate(Z[1]*0.8888888888888888 - Z[0]*0.4722222222222222) * W[1];
    // Sets the new, flowed SU3 matrix.
    m_updatedLattice[m_Index->getIndex(i,j,k,l)].U[mu].copy(exponentiate(m_S->getActionDerivative(lattice,W[2],i,j,k,l,mu)*m_epsilon*0.75 - Z[1]*0.8888888888888888 + Z[0]*0.4722222222222222)*W[2]);
    // HOW MUCH MEMORY THAT I WILL USE: (64**3*128*4*18*8)/1024/1024/1024*2 / (256/16)
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
     * Updates the lattice new values from m_updatedLattice.
     */
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        lattice[m_Index->getIndex(x,y,z,t)].U[mu].copy(m_updatedLattice[m_Index->getIndex(x,y,z,t)].U[mu]);
                    }
                }
            }
        }
    }
}
