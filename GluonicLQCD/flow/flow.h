#ifndef FLOW_H
#define FLOW_H

#include "math/lattice.h"
#include "math/flowexpfunctions.h"
#include "actions/action.h"

class Flow
{
private:
    // Flow update step
    double m_epsilon = 0.02;
    // Lattice constants
    std::vector<unsigned long int> m_N;
    unsigned long int m_subLatticeSize;
    // Temporary lattice to use when flowing
    Lattice<SU3> * m_tempLattice;
    Lattice<SU3> m_tempExpLattice;
    // Updates the lattice with the exponantiated lattice values
    inline Lattice<SU3> matrixExp(const Lattice<SU3> &lattice);
    inline Lattice<SU3> matrixExp(Lattice<SU3> &&lattice);
    // SU3 exponentiation function
    SU3Exp *m_SU3ExpFunc = nullptr;
    void setSU3ExpFunc();
    // Action pointer
    Action *m_S = nullptr;
public:
    Flow(Action *S);
    ~Flow();
    void flowField(Lattice<SU3> *lattice);
};

#endif // FLOW_H
