#include "mastersampler.h"

#include <cmath>
#include "parallelization/communicator.h"
#include "config/parameters.h"
#include "io/fieldio.h"

//using namespace LatticeOperations;

MasterSampler::MasterSampler(bool flow) : Correlator()
{
    // Sets up observable storage containers
    storeFlow(flow);

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
}

void MasterSampler::storeFlow(bool storeFlowObservable)
{
    m_storeFlowObservable = storeFlowObservable;
    if (m_storeFlowObservable) {
        // +1 as we are storing the initial value at t=0 as well.
        m_plaqObservable = new ObservableStorer(Parameters::getNFlows() + 1);
        m_topcObservable = new ObservableStorer(Parameters::getNFlows() + 1);
        m_energyObservable = new ObservableStorer(Parameters::getNFlows() + 1);
    } else {
        if (Parameters::getStoreThermalizationObservables()) {
            m_plaqObservable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
            m_topcObservable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
            m_energyObservable = new ObservableStorer(Parameters::getNCf() + Parameters::getNTherm() + 1);
        } else {
            m_plaqObservable = new ObservableStorer(Parameters::getNCf());
            m_topcObservable = new ObservableStorer(Parameters::getNCf());
            m_energyObservable = new ObservableStorer(Parameters::getNCf());
        }
    }
    m_plaqObservable->setNormalizeObservableByProcessor(true);
    m_plaqObservable->setObservableName("plaq");
    m_topcObservable->setObservableName("topc");
    m_energyObservable->setObservableName("energy");
}

void MasterSampler::writeFlowObservablesToFile(int iFlow)
{
    // Gathers and writes plaquette results to file
    m_plaqObservable->gatherResults();
    m_plaqObservable->writeFlowObservableToFile(iFlow);
    // Gathers and writes topological charge results to file
    m_topcObservable->gatherResults();
    m_topcObservable->writeFlowObservableToFile(iFlow);
    // Gathers and writes energy results to file
    m_energyObservable->gatherResults();
    m_energyObservable->writeFlowObservableToFile(iFlow);
}

void MasterSampler::writeObservableToFile(double acceptanceRatio)
{
    m_plaqObservable->writeObservableToFile(acceptanceRatio);
    m_topcObservable->writeObservableToFile(acceptanceRatio);
    m_energyObservable->writeObservableToFile(acceptanceRatio);
}

void MasterSampler::reset()
{
    /*
     * For resetting the flow observables between each flow.
     */
    m_plaqObservable->reset();
    m_topcObservable->reset();
    m_energyObservable->reset();

}

void MasterSampler::runStatistics()
{
    m_plaqObservable->runStatistics();
    m_topcObservable->runStatistics();
    m_energyObservable->runStatistics();
}

void MasterSampler::printHeader()
{
    if (!m_storeFlowObservable) {
        printf("%-*s %-*s %-*s",
               m_headerWidth,m_plaqObservable->getObservableName().c_str(),
               m_headerWidth,m_topcObservable->getObservableName().c_str(),
               m_headerWidth,m_energyObservable->getObservableName().c_str());
    } else {
        printf("\ni    t       %-*s %-*s %-*s",
               m_headerWidth,m_plaqObservable->getObservableName().c_str(),
               m_headerWidth,m_topcObservable->getObservableName().c_str(),
               m_headerWidth,m_energyObservable->getObservableName().c_str());
    }
}

void MasterSampler::printObservable(int iObs)
{
    if (!m_storeFlowObservable) {
        printf("%-*.8f %-*.8f %-*.8f",
               m_headerWidth,m_plaqObservable->getObservable(iObs),
               m_headerWidth,m_topcObservable->getObservable(iObs),
               m_headerWidth,m_energyObservable->getObservable(iObs));
    } else {
        double plaqObs = m_plaqObservable->getObservable(iObs); // TEMP TEMP TEMP!
        double topcObs = m_topcObservable->getObservable(iObs);
        double energyObs = m_energyObservable->getObservable(iObs);
        Parallel::Communicator::gatherDoubleResults(&plaqObs,1);
        Parallel::Communicator::gatherDoubleResults(&topcObs,1);
        Parallel::Communicator::gatherDoubleResults(&energyObs,1);
        Parallel::Communicator::setBarrier();
        if (Parallel::Communicator::getProcessRank() == 0) {
            printf("\n%-4d %-2.4f  %-*.15f %-*.15f %-*.15f",
                   iObs,
                   m_a*sqrt(8*Parameters::getFlowEpsilon()*iObs),
                   m_headerWidth,plaqObs/double(Parallel::Communicator::getNumProc()),
                   m_headerWidth,topcObs,
                   m_headerWidth,energyObs);
        }
    }
}

void MasterSampler::printStatistics()
{
    m_plaqObservable->printStatistics();
    m_topcObservable->printStatistics();
    m_energyObservable->printStatistics();
}

void MasterSampler::copyObservable(int iObs, std::vector<double> obs)
{
    m_plaqObservable[iObs] = obs[0];
    m_topcObservable[iObs] = obs[1];
    m_energyObservable[iObs] = obs[2];
}

std::vector<double> MasterSampler::getObservablesVector(int iObs)
{
    std::vector<double> obs(3);
    obs[0] = (*m_plaqObservable)[iObs];
    obs[1] = (*m_topcObservable)[iObs];
    obs[2] = (*m_energyObservable)[iObs];
    return obs;
}

void MasterSampler::calculate(Lattice<SU3> *lattice, int iObs)
{
    ///////////////////////////
    //// SYMMETRIC CLOVER /////
    ///////////////////////////
    m_topCharge = 0;
    m_energy = 0;
    m_plaquette = 0;
    mu = 0;

    for (int nu = 1; nu < 4; nu++)
    {
        // First clover. Definition from R Wohler 1985, more symmetric than other methods.
        // First leaf
        m_temp = lattice[mu];
        m_temp *= shift(lattice[nu],FORWARDS,mu);
        m_temp *= shift(lattice[mu],FORWARDS,nu).inv();
        m_temp *= lattice[nu].inv();
        m_clov1 = m_temp;

        // Adds plaquette leaf
        m_plaquette += sumRealTrace(m_clov1);

        // Retrieves beforehand in order to reduce number of communications by 2.
        m_U2Temp = shift(lattice[nu],BACKWARDS,nu);
        m_U3Temp = shift(lattice[mu],BACKWARDS,mu).inv();

        // Second leaf
        m_temp = lattice[mu];
        m_temp *= shift(shift(lattice[nu],FORWARDS,mu),BACKWARDS,nu).inv();
        m_temp *= shift(lattice[mu],BACKWARDS,nu).inv();
        m_temp *= m_U2Temp;
        m_clov1 -= m_temp;

        // Third leaf
        m_temp = m_U3Temp;
        m_temp *= shift(shift(lattice[nu],BACKWARDS,mu),BACKWARDS,nu).inv();
        m_temp *= shift(shift(lattice[mu],BACKWARDS,mu),BACKWARDS,nu);
        m_temp *= m_U2Temp;
        m_clov1 += m_temp;

        // Fourth leaf
        m_temp = m_U3Temp;
        m_temp *= shift(lattice[nu],BACKWARDS,mu);
        m_temp *= shift(shift(lattice[mu],FORWARDS,nu),BACKWARDS,mu);
        m_temp *= lattice[nu].inv();
        m_clov1 -= m_temp;

        int rho = nu % 3;
        rho++;
        int sigma = rho % 3;
        sigma++;

        // Second clover
        // First leaf
        m_temp = lattice[rho];
        m_temp *= shift(lattice[sigma],FORWARDS,rho);
        m_temp *= shift(lattice[rho],FORWARDS,sigma).inv();
        m_temp *= lattice[sigma].inv();
        m_clov2 = m_temp;

        // Gets lattice for temp use
        m_U2Temp = shift(lattice[sigma],BACKWARDS,sigma);
        m_U3Temp = shift(lattice[rho],BACKWARDS,rho).inv();

        // Adds another leaf to the plaquette
        m_plaquette += sumRealTrace(m_clov2);

        // Second leaf
        m_temp = lattice[rho];
        m_temp *= shift(shift(lattice[sigma],FORWARDS,rho),BACKWARDS,sigma).inv();
        m_temp *= shift(lattice[rho],BACKWARDS,sigma).inv();
        m_temp *= m_U2Temp;
        m_clov2 -= m_temp;

        // Third leaf
        m_temp = m_U3Temp;
        m_temp *= shift(shift(lattice[sigma],BACKWARDS,rho),BACKWARDS,sigma).inv();
        m_temp *= shift(shift(lattice[rho],BACKWARDS,rho),BACKWARDS,sigma);
        m_temp *= m_U2Temp;
        m_clov2 += m_temp;

        // Fourth leaf
        m_temp = m_U3Temp;
        m_temp *= shift(lattice[sigma],BACKWARDS,rho);
        m_temp *= shift(shift(lattice[rho],FORWARDS,sigma),BACKWARDS,rho);
        m_temp *= lattice[sigma].inv();
        m_clov2 -= m_temp;

        // Makes first clover anti hermitian and traceless
        m_temp = m_clov1.inv();
        m_clov1 -= m_temp;
        m_tempDiag = imagTrace(m_clov1)/3.0;
        m_clov1 = subtractImag(m_clov1,m_tempDiag);

        // Makes second clover anti hermitian and traceless
        m_temp = m_clov2.inv();
        m_clov2 -= m_temp;
        m_tempDiag = imagTrace(m_clov2)/3.0;
        m_clov2 = subtractImag(m_clov2,m_tempDiag);

        // Picks up the topological charge
        m_topCharge -= sumRealTraceMultiplication(m_clov1,m_clov2);

        // Picks up the action density
        m_energy += sumRealTraceMultiplication(m_clov1,m_clov1);
        m_energy += sumRealTraceMultiplication(m_clov2,m_clov2);
    }

    ///////////////////////////
    //////// PLAQUETTE ////////
    ///////////////////////////
    m_plaquette *= m_plaqMultiplicationFactor;
    (*m_plaqObservable)[iObs] = m_plaquette;
//    MPI_Allreduce(&m_plaquette,&m_plaquette,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    m_plaquette /= double(Parallel::Communicator::getNumProc());
//    if (Parallel::Communicator::getProcessRank() == 0) printf("\nPlaquette           = %20.16f",m_plaquette);

    ///////////////////////////
    //// TOPOLOGICAL CHARGE ///
    ///////////////////////////
    m_topCharge *= m_topcMultiplicationFactor;
    (*m_topcObservable)[iObs] = m_topCharge;
//    MPI_Allreduce(&m_topCharge,&m_topCharge,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    if (Parallel::Communicator::getProcessRank() == 0) printf("\nTopological charge  = %20.16f",m_topCharge);

    ///////////////////////////
    ///////// ENERGY //////////
    ///////////////////////////
    m_energy *= m_energyMultiplicationFactor;
    (*m_energyObservable)[iObs] = m_energy;
//    MPI_Allreduce(&m_energy,&m_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    if (Parallel::Communicator::getProcessRank() == 0) printf("\nEnergy              = %20.16f",m_energy);

}
