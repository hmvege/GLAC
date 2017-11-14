#include "energydensity.h"
#include "clover.h"
#include "math/functions.h"

EnergyDensity::EnergyDensity()
{
}

EnergyDensity::~EnergyDensity()
{

}

EnergyDensity::EnergyDensity(double a, int latticeSize) : Correlator()
{
    setLatticeSpacing(a);
    setLatticeSize(latticeSize);
    m_multiplicationFactor = (m_a*m_a*m_a*m_a)/(3*double(m_latticeSize));
}

void EnergyDensity::setLatticeSpacing(double a) // NEED TO DOUBLE CHECK THIS WITH ANDREA!
{
    m_multiplicationFactor = (a*a*a*a)/(3*double(m_latticeSize));
}

void EnergyDensity::calculate(SU3 *clovers, int i)
{
    m_actionDensity = 0;
    for (unsigned int i = 0; i < 12; i++)
    {
        m_actionDensity += traceSparseImagMultiplication(clovers[i],clovers[i]); // Might check this one with Andrea
    }
    return m_actionDensity;//*m_multiplicationFactor; // Correct or not?
}

void EnergyDensity::calculate(Links *lattice, int i)
{
    // When clover is not provided
    Clover Clov;
    Clov.setN(m_N);
    Clov.setLatticeSize(m_latticeSize);
    m_actionDensity = 0;
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    Clov.calculateClover(lattice,i,j,k,l);
                    for (unsigned int i = 0; i < 12; i++)
                    {
                        m_actionDensity += traceSparseImagMultiplication(Clov.m_clovers[i],Clov.m_clovers[i]);
                    }
                }
            }
        }
    }
    return m_actionDensity;//*m_multiplicationFactor; // Temporary off
}
