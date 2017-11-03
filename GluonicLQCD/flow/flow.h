#ifndef FLOW_H
#define FLOW_H

#include "links.h"

#include "actions/action.h"

#include "parallelization/indexorganiser.h"

class Flow
{
private:
    double m_epsilon = 0.01;
    double m_beta;
    unsigned int m_N[4];
    unsigned int m_subLatticeSize;
    Links * m_tempLattice;

    // Exponentiation variables for the Morningstar method
    complex f[3];
    complex h[3];
    double c0, c1, u, w, theta, xi0, c0max;
    SU3 f0; // Matrix for filling return matrix
    SU3 QSquared, QCubed, QQuartic;

    // Exponentiation variables for the Luscher method
//    SU3 U1,U2,U3,Y1,Y2,Y3,I,Y1Inv,Y2Inv,Y3Inv,E;
    SU3 U1,U2,U3,E,I;
    complex x1,x2,x3,div1,div2,div3;
    complex X1221X, X1331X, X2332X, sqrdFactor;

    // Exponentiation variables for the Taylor expansion method is just a reuse of QSquared and I

    // Parallelization
    int m_numprocs;
    int m_processRank;

    double cosu, cosw, sinu, sinw, uu, ww, sin2u, cos2u;

    // Functions for performing flow
    void runFlow(Links *lattice);
    IndexOrganiser *m_Index = nullptr;
    SU3 exponentiate(SU3 Q); // Morningstar method
    SU3 exponentiate2(SU3 Q); // Luscher method
    SU3 exponentiate3(SU3 Q); // Taylor expansion
    inline void updateLattice(Links *lattice);

    Action *m_S = nullptr;
public:
    Flow(unsigned int *N, double beta, int numprocs, int processRank);
    ~Flow();
    void flowGaugeField(int NFlows, Links *lattice);

    // Setters
    void setEpsilon(double epsilon) { m_epsilon = epsilon; }
    void setIndexHandler(IndexOrganiser *Index);
    void setAction(Action *S);
    // Getters
    double getEpsilon() { return m_epsilon; }
};

#endif // FLOW_H
