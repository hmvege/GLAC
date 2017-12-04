#include "clover.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include "math/functions.h"

Clover::Clover(bool storeFlowObservable) : Correlator(storeFlowObservable)
{
    /* The clover should be defined as:
     * index | mu | nu
     * 0       0    1
     * 1       0    2
     * 2       0    3
     * 3       1    2
     * 4       1    3
     * 5       2    3
     * 6       1    0
     * 7       2    0
     * 8       3    0
     * 9       2    1
     * 10      3    1
     * 11      3    2
     */
}

Clover::~Clover()
{

}

void Clover::calculateClover(Lattice<SU3> * lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    m_position[0] = i;
    m_position[1] = j;
    m_position[2] = k;
    m_position[3] = l;

    double tempDiag1, tempDiag2;
    int mu = 0;
    SU3 clov1, clov2;
    updateMuIndex(mu);
    for (int nu = 1; nu < 4; nu++)
    {
        updateNuIndex(nu);

        // First clover. Definition from R Wohler 1985, more symmetric than other methods.
        // First leaf
        U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
        U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
        clov1 = U1;

        // Adds plaquette leaf
        m_plaquettes[2*(nu-1)] = U1;

        // Retrieves beforehand in order to reduce number of communications by 2.
        U2Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
        U3Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,mu).inv();

        // Second leaf
        U2 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
        U2 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,mu).inv();
//        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
        U2 *= U2Temp;
        clov1 -= U2;


        // Third leaf
//        U3 = Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,mu).inv();
        U3 = U3Temp;
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,mu);
//        U3 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
        U3 *= U2Temp;
        clov1 += U3;

        // Fourth leaf
//        U4 = Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,mu).inv();
        U4 = U3Temp;
        U4 *= Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,nu);
        U4 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,nu,nuIndex,mu,muIndex,mu);
        U4 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
        clov1 -= U4;

        int rho = nu % 3;
        rho++;
        int sigma = rho % 3;
        sigma++;

        updateLorentzIndex(m_rhoIndex,rho);
        updateLorentzIndex(m_sigmaIndex,sigma);

        // Second clover
        // First leaf
        U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[rho];
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,rho,m_rhoIndex,sigma);
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,sigma,m_sigmaIndex,rho).inv();
        U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[sigma].inv();
        clov2 = U1;

        // Gets lattice for temp use
        U2Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,sigma,m_sigmaIndex,sigma);
        U3Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,rho,m_rhoIndex,rho).inv();

        m_plaquettes[2*nu - 1] = U1; // 2*(nu-1) + 1 = 2*nu - 2 + 1 = 2*nu - 1
        // Second leaf
        U2 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[rho];
        U2 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,rho,m_rhoIndex,sigma,m_sigmaIndex,sigma).inv();
        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,sigma,m_sigmaIndex,rho).inv();
//        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,sigma,m_sigmaIndex,sigma);
        U2 *= U2Temp;
        clov2 -= U2;

        // Third leaf
        U3 = U3Temp;
//        U3 = Parallel::Communicator::getNegativeLink(lattice,m_position,rho,m_rhoIndex,rho).inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,rho,m_rhoIndex,sigma,m_sigmaIndex,sigma).inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,rho,m_rhoIndex,sigma,m_sigmaIndex,rho);
//        U3 *= Parallel::Communicator::getNegativeLink(lattice,m_position,sigma,m_sigmaIndex,sigma);
        U3 *= U2Temp;
        clov2 += U3;

        // Fourth leaf
//        U4 = Parallel::Communicator::getNegativeLink(lattice,m_position,rho,m_rhoIndex,rho).inv();
        U4 = U3Temp;
        U4 *= Parallel::Communicator::getNegativeLink(lattice,m_position,rho,m_rhoIndex,sigma);
        U4 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,sigma,m_sigmaIndex,rho,m_rhoIndex,rho);
        U4 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[sigma].inv();
        clov2 -= U4;

        // Makes matrices anti hermitian and traceless
        U1 = clov1.inv();
        U2 = clov2.inv();
        clov1 -= U1;
        clov2 -= U2;
        tempDiag1 = (clov1[1] + clov1[9] + clov1[17])/3.0;
        tempDiag2 = (clov2[1] + clov2[9] + clov2[17])/3.0;
        for (int i = 1; i < 18; i+=8) {
            clov1[i] -= tempDiag1;
            clov2[i] -= tempDiag2;
        }

        m_clovers[2*(nu - 1)] = clov1;
        m_clovers[2*nu - 1] = clov2;
    }
}
