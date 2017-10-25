#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "links.h"
#include "parallelization/indexorganiser.h"
#include "parallelization/neighbours.h"
#include <vector>

class Correlator
{
protected:
    // Index length array
    unsigned int *m_N;
    double m_latticeSize;

    // For handling the shift-method in parallelization
    std::vector<int> indexes;
    IndexOrganiser *m_Index = nullptr;
public:
    Correlator();
    virtual ~Correlator();
    virtual double calculate(Links *lattice);
    virtual void setLatticeSize(int latticeSize);

    // Setters
    void setN(unsigned int *N);
    void initializeIndexHandler(IndexOrganiser *Index);
};


#endif // CORRELATOR_H
