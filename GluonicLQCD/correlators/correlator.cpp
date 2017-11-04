#include "correlator.h"
#include "parallelization/indexorganiser.h"
#include "parallelization/neighbours.h"

Correlator::Correlator()
{
    // Initiates the lattice dimensions
    m_N = new unsigned int[4];

    // Sets position vector to zero
    m_position = std::vector<int>(4,0);

    // Sets the lorentz indices to zero
    for (int i = 0; i < 4; i++) {
        muIndex[i] = 0;
        nuIndex[i] = 0;
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
    cout << "If you see this, something is wrong! Should not call correlator.cpp" << endl;
    return 1.0;
}

double Correlator::calculate()
{
    /*
     * Default correlator is not implemented.
     */
    cout << "If you see this, something is wrong! Should not call correlator.cpp" << endl;
    return 1.0;
}

void Correlator::setN(unsigned int *N) // MOVE INTO CONSTRUCTOR?
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}

void Correlator::initializeIndexHandler(IndexOrganiser *Index)
{
    m_Index = Index;
}
