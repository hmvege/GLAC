#include "correlator.h"
#include <iostream>

Correlator::Correlator(int N)
{
    m_N = N;
    m_latticeSize = N*N*N*N;
}

Correlator::~Correlator()
{

}

double Correlator::calculate(Links * lattice)
{
    /*
     * Default correlator is not implemented.
     */
    std::cout << "If you see this, something is wrong! Should not call correlator.cpp" << std::endl;
    return 1.0;
}
