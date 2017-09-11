#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"
#include "links.h"
#include "parallelization/indexorganiser.h"
#include <vector>

class Plaquette : public Correlator
{
private:
    int *muIndex;
    int *nuIndex;
    double multiplicationFactor;
    std::vector<int> indexes;
    IndexOrganiser *m_Index;
public:
//    Plaquette(int N, int N_T);
    Plaquette();
    ~Plaquette();
    double calculate(Links *lattice);
    void setLatticeSize(int latticeSize);
    void setIndexOrganiser(IndexOrganiser *Index) { m_Index = Index; }
};

#endif // PLAQUETTE_H
