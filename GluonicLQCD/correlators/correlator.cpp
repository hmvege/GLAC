#include "correlator.h"
#include "parallelization/indexorganiser.h"
#include "parallelization/neighbours.h"
#include <iostream>

Correlator::Correlator()
{
    m_N = new unsigned int[4];
    for (int i = 0; i < 4; i++) {
        m_N[i] = 0;
    }
    indexes = std::vector<int>(4,0);
}

void Correlator::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
}

Correlator::~Correlator()
{
    delete [] m_N;
}

double Correlator::calculate(Links * lattice)
{
    /*
     * Default correlator is not implemented.
     */
    std::cout << "If you see this, something is wrong! Should not call correlator.cpp" << std::endl;
    return 1.0;
}

void Correlator::setN(unsigned int *N)
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}

void Correlator::initializeIndexHandler(IndexOrganiser *Index)
{
    m_Index = Index;
}
