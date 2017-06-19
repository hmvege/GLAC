#ifndef ACTION_H
#define ACTION_H

#include "links.h"

class Action
{
protected:
    double a;
    double beta;
    int N;
public:
    Action();
    Action(int latticeSize, double new_a, double new_beta);
    virtual ~Action() {}
    virtual double getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu);
    // Setters
    void setLatticeSpacing(double new_a) { a = new_a; }
    void setNLatticePoints(int new_N) { N = new_N; }
    // Getters
    double getLatticeSpacing() { return a; }
    int getNLatticePoints() { return N; }
};

inline int lorentzUnit(int mu) {
    /*
     * Unit vector for lorentz indexes
     */
    if (mu==0) {
        return N*N*N;
    } else if (mu==1) {
        return N*N;
    } else if (mu==2) {
        return N;
    } else {
        return 1;
    }
}

#endif // ACTION_H
