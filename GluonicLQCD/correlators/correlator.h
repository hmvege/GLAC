#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "links.h"
#include "parallelization/indexorganiser.h"
#include <vector>

class Correlator
{
protected:
    // Index length array
    int *m_N;
    double m_latticeSize;
    // For handling the shift-method in parallelization
    std::vector<int> indexes;
    IndexOrganiser *m_Index;
public:
//    Correlator(int N, int N_T);
    Correlator();
    virtual ~Correlator();
    virtual double calculate(Links *lattice);
    virtual void setLatticeSize(int latticeSize);
    // Setters
    void setN(int *N);
    void setIndexOrganiser(IndexOrganiser *Index) { m_Index = Index; }
};

#endif // CORRELATOR_H
