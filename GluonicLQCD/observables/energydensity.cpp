#include "energydensity.h"

EnergyDensity::EnergyDensity(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    m_observable->setObservableName(m_observableNameCompact);
    m_observable->setNormalizeObservableByProcessor(false);

    // Allocates disc space for the lattices to be used for calculations.
    m_clov1.allocate(m_N);
    m_clov2.allocate(m_N);
    m_U2Temp.allocate(m_N);
    m_U3Temp.allocate(m_N);
    m_temp.allocate(m_N);

    // Sets observable multiplication factor
    m_multiplicationFactor = 1.0/double(m_latticeSize);
}

EnergyDensity::~EnergyDensity()
{

}

void EnergyDensity::calculate(Lattice<SU3> *lattice, unsigned int iObs)
{
    m_energyDensity = 0;
    mu = 0;

    for (int nu = 1; nu < 4; nu++)
    {
        // First clover. Definition from R Wohler 1985, more symmetric than other methods.
        // First leaf
        m_temp = lattice[mu];
        m_temp *= shift(lattice[nu],FORWARDS,mu);
        m_temp *= inv(shift(lattice[mu],FORWARDS,nu));
        m_temp *= inv(lattice[nu]);
        m_clov1 = m_temp;

        // Retrieves beforehand in order to reduce number of communications by 2.
        m_U2Temp = shift(lattice[nu],BACKWARDS,nu);
        m_U3Temp = inv(shift(lattice[mu],BACKWARDS,mu));

        // Second leaf
        m_temp = lattice[mu];
        m_temp *= inv(shift(shift(lattice[nu],FORWARDS,mu),BACKWARDS,nu));
        m_temp *= inv(shift(lattice[mu],BACKWARDS,nu));
        m_temp *= m_U2Temp;
        m_clov1 -= m_temp;

        // Third leaf
        m_temp = m_U3Temp;
        m_temp *= inv(shift(shift(lattice[nu],BACKWARDS,mu),BACKWARDS,nu));
        m_temp *= shift(shift(lattice[mu],BACKWARDS,mu),BACKWARDS,nu);
        m_temp *= m_U2Temp;
        m_clov1 += m_temp;

        // Fourth leaf
        m_temp = m_U3Temp;
        m_temp *= shift(lattice[nu],BACKWARDS,mu);
        m_temp *= shift(shift(lattice[mu],FORWARDS,nu),BACKWARDS,mu);
        m_temp *= inv(lattice[nu]);
        m_clov1 -= m_temp;

        // Updating indexes for a second clover
        rho = nu % 3;
        rho++;
        sigma = rho % 3;
        sigma++;

        // Second clover
        // First leaf
        m_temp = lattice[rho];
        m_temp *= shift(lattice[sigma],FORWARDS,rho);
        m_temp *= inv(shift(lattice[rho],FORWARDS,sigma));
        m_temp *= inv(lattice[sigma]);
        m_clov2 = m_temp;

        // Gets lattice for temp use
        m_U2Temp = shift(lattice[sigma],BACKWARDS,sigma);
        m_U3Temp = inv(shift(lattice[rho],BACKWARDS,rho));

        // Second leaf
        m_temp = lattice[rho];
        m_temp *= inv(shift(shift(lattice[sigma],FORWARDS,rho),BACKWARDS,sigma));
        m_temp *= inv(shift(lattice[rho],BACKWARDS,sigma));
        m_temp *= m_U2Temp;
        m_clov2 -= m_temp;

        // Third leaf
        m_temp = m_U3Temp;
        m_temp *= inv(shift(shift(lattice[sigma],BACKWARDS,rho),BACKWARDS,sigma));
        m_temp *= shift(shift(lattice[rho],BACKWARDS,rho),BACKWARDS,sigma);
        m_temp *= m_U2Temp;
        m_clov2 += m_temp;

        // Fourth leaf
        m_temp = m_U3Temp;
        m_temp *= shift(lattice[sigma],BACKWARDS,rho);
        m_temp *= shift(shift(lattice[rho],FORWARDS,sigma),BACKWARDS,rho);
        m_temp *= inv(lattice[sigma]);
        m_clov2 -= m_temp;

        // Makes first clover anti hermitian and traceless
        m_temp = inv(m_clov1);
        m_clov1 -= m_temp;
        m_tempDiag = imagTrace(m_clov1)/3.0;
        m_clov1 = subtractImag(m_clov1,m_tempDiag);

        // Makes second clover anti hermitian and traceless
        m_temp = inv(m_clov2);
        m_clov2 -= m_temp;
        m_tempDiag = imagTrace(m_clov2)/3.0;
        m_clov2 = subtractImag(m_clov2,m_tempDiag);

        // Picks up the energy density
        m_energyDensity += sumRealTraceMultiplication(m_clov1,m_clov1);
        m_energyDensity += sumRealTraceMultiplication(m_clov2,m_clov2);
    }
    (*m_observable)[iObs] = m_energyDensity*m_multiplicationFactor;
}

void EnergyDensity::runStatistics()
{
    m_observable->gatherResults();
    m_observable->runStatistics();
}

void EnergyDensity::printStatistics()
{
    m_observable->printStatistics();
}
