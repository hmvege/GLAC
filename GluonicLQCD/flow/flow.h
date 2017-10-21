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

    // Flow variables
    complex f[3];
    complex h[3];
    complex c0, c1, u, w, theta, xi0;
    SU3 W[3];
    SU3 T[8]; // Gell-Mann matrices made anti-hermitian and traceless by multiplying with (-i). Might not need this,a s the Gell-Mann matrices sets so many numbers to zero?
    SU3 I; // Identity
    SU3 QSquared, QCubed;

    // Functions for performing flow
    void runFlow(Links *lattice);
    void smearLink(SU3 V);
    IndexOrganiser *m_Index = nullptr;
    SU3 derivative(Links * lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    SU3 exponentiate(SU3 X);

    Action *m_S = nullptr;
public:
    Flow(unsigned int *N, double beta);
    void flowGaugeField(int NFlows, Links *lattice);

    // Setters
    void setEpsilon(double epsilon) { m_epsilon = epsilon; }
    void setIndexHandler(IndexOrganiser *Index);
    void setAction(Action *S) { m_S = S; }
    // Getters
    double getEpsilon() { return m_epsilon; }
};

#endif // FLOW_H
