#include "correlator.h"

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
    m_observable = new ObservableStorer(NSize);
}

void Correlator::setLatticeSize(int latticeSize)
{
    m_latticeSize = double(latticeSize);
}

Correlator::~Correlator()
{
    delete [] m_N;
}

void Correlator::calculate(Links * lattice, int i)
{
    /*
     * Default correlator is not implemented.
     */
    lattice[i].U[0].zeros(); // TEMP
    printf("\nIf you see this, something is wrong! Should not call correlator.cpp");
}

void Correlator::calculate(SU3 *U, int i)
{
    /*
     * Default correlator is not implemented when only given a SU3 matrix.
     */
    U[i % 4].zeros(); // TEMP
    printf("\nIf you see this, something is wrong! Should not call correlator.cpp");
}

void Correlator::setN(unsigned int *N) // MOVE INTO CONSTRUCTOR?
{
    for (int i = 0; i < 4; i++) {
        m_N[i] = N[i];
    }
}

double Correlator::getObservable(int i)
{
    printf("\nNot implemented for base class!");
    return i;
}

void Correlator::writeStatisticsToFile()
{
    printf("\nFunction for writing statistics to file not implemented for base correlator class!");
}
