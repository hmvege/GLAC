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
    Links * m_updatedLattice;

    // Flow variables
    complex f[3];
    complex h[3];
    double c0, c1, u, w, theta, xi0, c0max;
    SU3 W[3], Z[2];
    SU3 I; // Matrix for filling return matrix
    SU3 QSquared, QCubed;

    double temp, cosu, cosw, sinu, sinw, uu, ww, ww3uu, sin2u, cos2u;

    // Functions for performing flow
    void runFlow(Links *lattice);
    void smearLink(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu);
    IndexOrganiser *m_Index = nullptr;
    SU3 exponentiate(SU3 Q);
    inline void updateLattice(Links *lattice);

    Action *m_S = nullptr;
public:
    Flow(unsigned int *N, double beta);
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
