#ifndef ACTION_H
#define ACTION_H

#include "links.h"
#include "functions.h"

class Action
{
protected:
    double a;
    double beta;
    int N;
    // Lorentz indices arrays
    int *muIndex;
    int *nuIndex;
public:
    Action();
    Action(int latticeSize, double new_a, double new_beta);
    virtual ~Action();
    virtual double getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu);
    // Setters
    void setLatticeSpacing(double new_a) { a = new_a; }
    void setNLatticePoints(int new_N) { N = new_N; }
    // Getters
    double getLatticeSpacing() { return a; }
    int getNLatticePoints() { return N; }
};

inline int stapleIndex(int i, int j, int k, int l, int N)
{
    /*
     * Unit vector for lorentz indexes
     */
    return index((i+N) % N, (j+N) % N, (k+N) % N, (l+N) % N, N);;
}

inline void lorentzIndex(int mu, int *lorentzIndices)
{
    /*
     * Fills the lorentz array with correct indices
     */
    for (int i = 0; i < 4; i++)
    {
        if (mu==i)
        {
            lorentzIndices[i] = 1;
        }
        else
        {
            lorentzIndices[i] = 0;
        }
    }
}


#endif // ACTION_H
