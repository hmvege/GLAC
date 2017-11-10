#include "clover.h"

Clover::Clover() : Correlator()
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
    m_overCounter = 0; // Dirty method of ensuring cloverIndex returns right value.
    for (int mu = 0; mu < 4; mu++)
    {
        updateMuIndex(mu);
        for (int nu = mu; nu < 4; nu++) // ONLY NEED TO STORE HALF OF THESE; SINCE THEY ARE SYMMETRIC AND CAN BE TAKEN INVERSE OF!
        {
            if (nu==mu) {
                m_overCounter++;
                continue;
            }
            updateNuIndex(nu);

            // First leaf
            U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
            U1 *= Parallel::Index::getPositiveLink(lattice,m_position,mu,muIndex,nu);
            U1 *= Parallel::Index::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
            U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();

            // Second leaf
            U2 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu];
            U2 *= Parallel::Index::getNeighboursNeighbourLink(lattice,m_position,nu,nuIndex,mu,muIndex,mu).inv();
            U2 *= Parallel::Index::getNegativeLink(lattice,m_position,mu,muIndex,nu).inv();
            U2Temp = Parallel::Index::getNegativeLink(lattice,m_position,mu,muIndex,mu);
            U2 *= U2Temp;

            // Third leaf
            U3 = U2Temp.inv();
            U3 *= Parallel::Index::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
            U3 *= Parallel::Index::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,mu);
            U3Temp = Parallel::Index::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
            U3 *= U3Temp;

            // Fourth leaf
            U4 = U3Temp.inv();
            U4 *= Parallel::Index::getNegativeLink(lattice,m_position,nu,nuIndex,mu);
            U4 *= Parallel::Index::getNeighboursNeighbourLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu);
            U4 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu].inv();

            // Gets the plaquette leaf
//            printf("mu=%2d nu=%2d plaq=%2d clov=%2d clov_inv=%2d\n",mu,nu,3*mu + nu - m_overCounter,cloverIndex(mu,nu-m_overCounter),3*nu+mu);
            m_plaquettes[3*mu + nu - m_overCounter - mu/2] = U1;
            // Sums the leafs, takes imaginary part(sets real values to zero) and multiply by 0.25.
//            m_clovers[cloverIndex(mu,nu-m_overCounter)] = (U1 + U2 + U3 + U4).getIm()*0.25;
//            m_clovers[3*nu + mu] = m_clovers[cloverIndex(mu,nu-m_overCounter)].inv();

            SU3 A;
            A = (U1 + U2 + U3 + U4);
            m_clovers[cloverIndex(mu,nu-m_overCounter)] = (A - A.inv()) * (1/8.0); // Using the old luscher definition

            m_clovers[3*nu + mu] = m_clovers[cloverIndex(mu,nu-m_overCounter)].inv();
//            m_clovers[cloverIndex(mu,nu-m_overCounter)].printMachine();
//            exit(1);
        }
    }
}
