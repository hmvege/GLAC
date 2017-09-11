#include "correlator.h"
#include <iostream>


Correlator::Correlator()
{
    m_N = new int[4];
    for (int i = 0; i < 4; i++) {
        m_N[i] = 0;
    }
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

void Correlator::setN(int *N)
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}
