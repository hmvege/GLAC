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

void Clover::calculateClover(Links *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    m_position[0] = i;
    m_position[1] = j;
    m_position[2] = k;
    m_position[3] = l;
//    // OLD METHOD
//    m_overCounter = 0; // Dirty method of ensuring cloverIndex returns right value.
//    for (int mu = 0; mu < 4; mu++)
//    {
//        updateMuIndex(mu);
//        for (int nu = mu; nu < 4; nu++) // ONLY NEED TO STORE HALF OF THESE; SINCE THEY ARE SYMMETRIC AND CAN BE TAKEN INVERSE OF!
//        {
//            if (nu==mu) {
//                m_overCounter++;
//                continue;
//            }
//            updateNuIndex(nu);

//            // First leaf
//            U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
//            U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
//            U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
//            U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();

////            // Second leaf
////            U2 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu];
////            U2 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,nu,nuIndex,mu,muIndex,mu).inv();
////            U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,nu).inv();
////            U2Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,mu);
////            U2 *= U2Temp;

////            // Third leaf
////            U3 = U2Temp.inv();
////            U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
////            U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,mu);
////            U3Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
////            U3 *= U3Temp;

////            // Fourth leaf
////            U4 = U3Temp.inv();
////            U4 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,mu);
////            U4 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu);
////            U4 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu].inv();

//            // Gets the plaquette leaf
////            printf("mu=%2d nu=%2d plaq=%2d clov=%2d clov_inv=%2d\n",mu,nu,3*mu + nu - m_overCounter,cloverIndex(mu,nu-m_overCounter),3*nu+mu);
//            m_plaquettes[3*mu + nu - m_overCounter - mu/2] = U1;
//            // Sums the leafs, takes imaginary part(sets real values to zero) and multiply by 0.25.
////            m_clovers[cloverIndex(mu,nu-m_overCounter)] = (U1 + U2 + U3 + U4).getIm()*0.25;
////            m_clovers[3*nu + mu] = m_clovers[cloverIndex(mu,nu-m_overCounter)].inv();

////            SU3 A;
////            A = (U1 + U2 + U3 + U4);
////            m_clovers[cloverIndex(mu,nu-m_overCounter)] = (A - A.inv()) * (1/8.0); // Using the old luscher definition
//            m_clovers[cloverIndex(mu,nu-m_overCounter)] = (U1 - U1.inv()) * (1/8.0);
//            m_clovers[3*nu + mu] = m_clovers[cloverIndex(mu,nu-m_overCounter)].inv();
////            m_clovers[cloverIndex(mu,nu-m_overCounter)].printMachine();
////            exit(1);
//        }
//    }

    int rhoIndex[4];
    int sigmaIndex[4];
    // CHROMA METHOD
    int mu = 0;
    SU3 clov1, clov2,I;
    I.identity();
    updateMuIndex(mu);
    m_overCounter = 0; // Dirty method of ensuring cloverIndex returns right value.
    for (int nu = 1; nu < 4; nu++)
    {
        updateNuIndex(nu);

        // First clover
        // First leaf
        U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
        U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
        clov1 = U1;

        m_plaquettes[2*(nu-1)] = U1;
//        printf("\n nu=%d      2*(nu-1) = %d", nu, 2*(nu-1));
        // Second leaf
        U2 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu];
        U2 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,nu,nuIndex,mu,muIndex,mu).inv();
        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,nu).inv();
        U2Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,mu);
        U2 *= U2Temp;
        clov1 -= U2;

        // Third leaf
        U3 = U2Temp.inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,mu);
        U3Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
        U3 *= U3Temp;
        clov1 += U3;

        // Fourth leaf
        U4 = U3Temp.inv();
        U4 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,mu);
        U4 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu);
        U4 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu].inv();
        clov1 -= U4;

        int rho = nu % 3;
        rho++;
        int sigma = rho % 3;
        sigma++;

        updateLorentzIndex(rhoIndex,rho);
        updateLorentzIndex(sigmaIndex,sigma);

        // Second clover
        // First leaf
        U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[rho];
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,rho,rhoIndex,sigma);
        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,sigma,sigmaIndex,rho).inv();
        U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[sigma].inv();
        clov2 = U1;

        m_plaquettes[2*nu - 1] = U1; // 2*(nu-1) + 1 = 2*nu - 2 + 1 = 2*nu - 1
//        printf("\n nu=%d      2*nu - 1 = %d", nu, 2*nu - 1);

        // Second leaf
        U2 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[sigma];
        U2 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,sigma,sigmaIndex,rho,rhoIndex,rho).inv();
        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,rho,rhoIndex,sigma).inv();
        U2Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,rho,rhoIndex,rho);
        U2 *= U2Temp;
        clov2 -= U2;

        // Third leaf
        U3 = U2Temp.inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,rho,rhoIndex,sigma,sigmaIndex,sigma).inv();
        U3 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,rho,rhoIndex,sigma,sigmaIndex,rho);
        U3Temp = Parallel::Communicator::getNegativeLink(lattice,m_position,sigma,sigmaIndex,sigma);
        U3 *= U3Temp;
        clov2 += U3;

        // Fourth leaf
        U4 = U3Temp.inv();
        U4 *= Parallel::Communicator::getNegativeLink(lattice,m_position,sigma,sigmaIndex,rho);
        U4 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,rho,rhoIndex,sigma,sigmaIndex,sigma);
        U4 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[rho].inv();
        clov2 -= U4;

        clov1 = (clov1 - clov1.inv()) - I*clov1.trace()/3.0;
        clov2 = (clov2 - clov2.inv()) - I*clov2.trace()/3.0;

        m_clovers[2*(nu - 1)] = clov1;
        m_clovers[2*nu - 1] = clov2;
//        clov1.printMachine();
//        clov2.printMachine();
    }
//    Parallel::Communicator::MPIExit("\nExiting at clover.cpp");
}
