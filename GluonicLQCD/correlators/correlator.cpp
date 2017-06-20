#include "correlator.h"

Correlator::Correlator()
{

}

Correlator::Correlator(int latticeSize)
{
    m_latticeSize = latticeSize;
}

Correlator::~Correlator()
{

}

double Correlator::calculate(Links * lattice)
{
    return 1.0;
}
