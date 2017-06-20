#include "correlator.h"

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
    return 1.0;
}
