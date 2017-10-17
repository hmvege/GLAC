#ifndef FLOW_H
#define FLOW_H

#include "links.h"

#include "parallelization/indexorganiser.h"

class Flow
{
private:
    double m_epsilon = 0.01;
    int m_N[4];
    // Functions for performing flow
    void runFlow(Links *lattice);
    void smearLink();
    IndexOrganiser *m_Index = nullptr;
public:
    Flow();
    void flowGaugeField(int NFlows, Links *lattice);

    // Setters
    void setEpsilon(double epsilon) { m_epsilon = epsilon; }
    // Getters
    double getEpsilon() { return m_epsilon; }
};

#endif // FLOW_H
