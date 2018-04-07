#include "mastersamplertopcxyz.h"

#include <cmath>
#include "parallelization/communicator.h"
#include "config/parameters.h"
#include "io/observablesio.h"

MasterSamplerTopcXYZ::MasterSamplerTopcXYZ(bool flow) : Correlator()
{
    // Sets up observable storage containers
    initializeObservableStorer(flow);

    // Sets up multiplication factors
    m_plaqMultiplicationFactor = 1.0/(18.0*double(m_latticeSize));
    m_topcMultiplicationFactor = 1.0/(16*16*M_PI*M_PI);
    m_energyMultiplicationFactor = 1.0/double(Parameters::getLatticeSize()); // Cant divide by size if we are not normalizing as well

    // Allocates memory to the helper variables
    m_clov1.allocate(m_N);
    m_clov2.allocate(m_N);
    m_U2Temp.allocate(m_N);
    m_U3Temp.allocate(m_N);
    m_temp.allocate(m_N);

    // Allocates temporary vector for retrieves results from the lattice method
    m_tempTopct.resize(m_N[3]);

    // Allocates temporary array for gathering results into a single time array
    m_topctGatherVector.resize(Parameters::getNTemporal() * (Parameters::getNFlows() + 1));
    for (unsigned int iFlow = 0; iFlow < Parameters::getNFlows() + 1; iFlow++) {
        for (unsigned int it = 0; it < Parameters::getNTemporal(); it++) {
            m_topctGatherVector[iFlow*Parameters::getNTemporal() + it] = 0;
        }
    }
}

MasterSamplerTopcXYZ::~MasterSamplerTopcXYZ()
{
    delete m_plaqObservable;
    delete m_topcObservable;
    delete m_energyObservable;
    delete m_topctObservable;
}

void MasterSamplerTopcXYZ::initializeObservableStorer(bool storeFlowObservable)
{
    m_storeFlowObservable = storeFlowObservable;
    if (m_storeFlowObservable) {
        // +1 as we are storing the initial value at t=0 as well.
        m_plaqObservable = new ObservableStorer(Parameters::getNFlows() + 1);
        m_topcObservable = new ObservableStorer(Parameters::getNFlows() + 1);
        m_energyObservable = new ObservableStorer(Parameters::getNFlows() + 1);
        m_topctObservable = new ObservableStorer((Parameters::getNFlows() + 1) * m_N[3]);
    } else {
        if (Parameters::getStoreThermalizationObservables()) {
            m_plaqObservable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
            m_topcObservable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
            m_energyObservable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
            m_topctObservable = new ObservableStorer((Parameters::getNCf() + Parameters::getNTherm() + 1) * m_N[3]);
        } else {
            m_plaqObservable = new ObservableStorer(Parameters::getNCf());
            m_topcObservable = new ObservableStorer(Parameters::getNCf());
            m_energyObservable = new ObservableStorer(Parameters::getNCf());
            m_topctObservable = new ObservableStorer(Parameters::getNCf() * m_N[3]);
        }
    }
    m_plaqObservable->setNormalizeObservableByProcessor(true);
    m_plaqObservable->setObservableName("plaq");
    m_topcObservable->setObservableName("topc");
    m_energyObservable->setObservableName("energy");
    m_topctObservable->setObservableName("topct");
}

void MasterSamplerTopcXYZ::writeFlowObservablesToFile(unsigned int configNumber)
{
    // Gathers and writes plaquette results to file
    m_plaqObservable->gatherResults();
    m_plaqObservable->writeFlowObservableToFile(configNumber);

    // Gathers and writes topological charge results to file
    m_topcObservable->gatherResults();
    m_topcObservable->writeFlowObservableToFile(configNumber);

    // Gathers and writes energy results to file
    m_energyObservable->gatherResults();
    m_energyObservable->writeFlowObservableToFile(configNumber);

    // Flattens topct to t direction and gathers topct results into a single array
    Parallel::Communicator::reduceToTemporalDimension(m_topctGatherVector, m_topctObservable->getObservableArray());

    // Gathers and writes the Euclidean time array
    IO::writeMatrixToFile(m_topctGatherVector,m_topctObservable->getObservableName(),configNumber,Parameters::getNTemporal());
}

void MasterSamplerTopcXYZ::writeObservableToFile(double acceptanceRatio)
{
    // Writes plaquette results to file
    m_plaqObservable->writeObservableToFile(acceptanceRatio);

    // Writes topological charge results to file
    m_topcObservable->writeObservableToFile(acceptanceRatio);

    // Writes energy results to file
    m_energyObservable->writeObservableToFile(acceptanceRatio);
}

void MasterSamplerTopcXYZ::reset()
{
    /*
     * For resetting the flow observables between each flow.
     */
    m_plaqObservable->reset();
    m_topcObservable->reset();
    m_energyObservable->reset();
    m_topctObservable->reset();
    for (unsigned int i = 0; i < Parameters::getNTemporal() * (Parameters::getNFlows() + 1); i++) {
        m_topctGatherVector[i] = 0;
    }

//    m_clov1.zeros();
//    m_clov2.zeros();
//    m_U2Temp.zeros();
//    m_U3Temp.zeros();
//    m_temp.zeros();
}

void MasterSamplerTopcXYZ::runStatistics()
{
    m_plaqObservable->gatherResults();
    m_plaqObservable->runStatistics();

    m_topcObservable->gatherResults();
    m_topcObservable->runStatistics();

    m_energyObservable->gatherResults();
    m_energyObservable->runStatistics();
}

void MasterSamplerTopcXYZ::printHeader()
{
    if (!m_storeFlowObservable) {
        printf("%-*s %-*s %-*s %-*s",
               m_headerWidth,m_plaqObservable->getObservableName().c_str(),
               m_headerWidth,m_topcObservable->getObservableName().c_str(),
               m_headerWidth,m_energyObservable->getObservableName().c_str(),
               m_headerWidth,m_topctObservable->getObservableName().c_str());
    } else {
        printf("\ni    t      %-*s %-*s %-*s %-*s",
               m_headerWidth,m_plaqObservable->getObservableName().c_str(),
               m_headerWidth,m_topcObservable->getObservableName().c_str(),
               m_headerWidth,m_energyObservable->getObservableName().c_str(),
               m_headerWidth,m_topctObservable->getObservableName().c_str());
    }
}

void MasterSamplerTopcXYZ::printObservable(unsigned int iObs)
{
    if (!m_storeFlowObservable) {
        double topctObs = 0;
        for (unsigned int it = 0; it < m_N[3]; it++) {
            topctObs += (*m_topctObservable)[iObs*m_N[3] + it];
        }

        double topcObs = 0;
        topcObs = m_topcObservable->getObservable(iObs);

        if (Parameters::getNFlows() != 0) {
            Parallel::Communicator::gatherDoubleResults(&topcObs,1);
        }
        Parallel::Communicator::gatherDoubleResults(&topctObs,1);

        if (Parallel::Communicator::getProcessRank() == 0) {
            printf("%-*.8f %-*.8f %-*.8f %-*.8f",
                   m_headerWidth,m_plaqObservable->getObservable(iObs),
                   m_headerWidth,topcObs,
                   m_headerWidth,m_energyObservable->getObservable(iObs),
                   m_headerWidth,topctObs);
        }
    } else {
        double plaqObs = m_plaqObservable->getObservable(iObs);
        double topcObs = m_topcObservable->getObservable(iObs);
        double energyObs = m_energyObservable->getObservable(iObs);
        double topctObs = 0;
        for (unsigned int it = 0; it < m_N[3]; it++) {
            topctObs += (*m_topctObservable)[iObs*m_N[3] + it];
        }

        Parallel::Communicator::gatherDoubleResults(&plaqObs,1);
        Parallel::Communicator::gatherDoubleResults(&topcObs,1);
        Parallel::Communicator::gatherDoubleResults(&energyObs,1);
        Parallel::Communicator::gatherDoubleResults(&topctObs,1);

        if (Parallel::Communicator::getProcessRank() == 0) {
            printf("\n%-4d %-3.3f  %-*.15f %-*.15f %-*.15f %-*.15f",
                   iObs,
                   double(iObs)*Parameters::getFlowEpsilon(),
                   m_headerWidth,plaqObs/double(Parallel::Communicator::getNumProc()),
                   m_headerWidth,topcObs,
                   m_headerWidth,energyObs,
                   m_headerWidth,topctObs);
        }
    }
}

void MasterSamplerTopcXYZ::printStatistics()
{
    m_plaqObservable->printStatistics();
    m_topcObservable->printStatistics();
    m_energyObservable->printStatistics();
}

void MasterSamplerTopcXYZ::copyObservable(unsigned int iObs, std::vector<double> obs)
{
    (*m_plaqObservable)[iObs] = obs[0];
    (*m_topcObservable)[iObs] = obs[1];
    (*m_energyObservable)[iObs] = obs[2];
    for (unsigned long int it = 0; it < m_N[3]; it++) {
        (*m_topctObservable)[iObs*m_N[3] + it] = obs[3 + it];
    }
}

std::vector<double> MasterSamplerTopcXYZ::getObservablesVector(unsigned int iObs)
{
    std::vector<double> obs(3 + m_N[3]);
    obs[0] = (*m_plaqObservable)[iObs];
    obs[1] = (*m_topcObservable)[iObs];
    obs[2] = (*m_energyObservable)[iObs];
    for (unsigned long int it = 0; it < m_N[3]; it++) {
        obs[3 + it] = (*m_topctObservable)[iObs*m_N[3] + it];
    }

    return obs;
}

void MasterSamplerTopcXYZ::calculate(Lattice<SU3> *lattice, unsigned int iObs)
{
    ///////////////////////////
    //// SYMMETRIC CLOVER /////
    ///////////////////////////
    m_topCharge = 0;
    m_energy = 0;
    m_plaquette = 0;
    for (unsigned int it = 0; it < m_N[3]; it++) {
        m_tempTopct[it] = 0;
    }
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

        // Adds plaquette leaf
        m_plaquette += sumRealTrace(m_clov1);

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

        // Adds another leaf to the plaquette
        m_plaquette += sumRealTrace(m_clov2);

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

        // Sums take the real trace multiplication and sums into a temporary holder
        m_tempTopct = sumXYZ(realTraceMultiplication(m_clov1,m_clov2));

        // Loops over time dimension
        for (unsigned long int it = 0; it < m_N[3]; it++) {
            // Sums the topological charge values at the xyz axis into a observable holder
            (*m_topctObservable)[iObs*m_N[3] + it] -= m_tempTopct[it];

            // Sums the topological charge, negative sign already taken care of
            m_topCharge -= m_tempTopct[it];
        }

        // Picks up the action density
        m_energy += sumRealTraceMultiplication(m_clov1,m_clov1);
        m_energy += sumRealTraceMultiplication(m_clov2,m_clov2);
    }


    ///////////////////////////
    //////// PLAQUETTE ////////
    ///////////////////////////
    m_plaquette *= m_plaqMultiplicationFactor;
    (*m_plaqObservable)[iObs] = m_plaquette;

    ///////////////////////////
    //// TOPOLOGICAL CHARGE ///
    ///////////////////////////
    m_topCharge *= m_topcMultiplicationFactor;
    (*m_topcObservable)[iObs] = m_topCharge;

    ///////////////////////////
    ///////// ENERGY //////////
    ///////////////////////////
    m_energy *= m_energyMultiplicationFactor;
    (*m_energyObservable)[iObs] = m_energy;

    ///////////////////////////
    ///// TOP. CHARGE XYZ /////
    ///////////////////////////
    for (unsigned int it = 0; it < m_N[3]; it++) {
        (*m_topctObservable)[iObs*m_N[3] + it] *= m_topcMultiplicationFactor;
    }
}
