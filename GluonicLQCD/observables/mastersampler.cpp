#include "mastersampler.h"

#include "parallelization/communicator.h"
#include "config/parameters.h"
#include "io/fieldio.h"

//using namespace LatticeOperations;

MasterSampler::MasterSampler()
{
    m_multiplicationFactor = 18.0*double(Parameters::getLatticeSize());

}

void MasterSampler::calculate()
{
    // Initializes lattice
    Lattice <SU3> lattice[4];
    unsigned int N[4] = {4, 8, 8, 16};
    Parameters::setSubLatticePreset(true);
    Parameters::setN(N);
    Parallel::Index::setN(N);
    Parallel::Communicator::setN(N);
    Parallel::Communicator::initializeSubLattice();
    m_multiplicationFactor = 18.0*double(Parameters::getSubLatticeSize());
    IO::FieldIO::init();
    std::vector<int> dim = {int(N[0]),int(N[1]),int(N[2]),int(N[3])};
    printf("\n%d %d %d %d \n",int(N[0]),int(N[1]),int(N[2]),int(N[3]));
    for (int mu = 0; mu < 4; mu++) {
        lattice[mu].allocate(dim);
        lattice[mu].identity();
    }

    // Loads configuration into lattice
    std::string fname = "LatticeOperationsTestConfig_beta6.000000_spatial8_temporal16_threads4_config0.bin"; // 0.593424
    IO::FieldIO::loadLatticeFieldConfiguration(fname,lattice);
    for (int mu = 0; mu < 4; mu++) {
        for (int iSite = 0; iSite < lattice[mu].m_latticeSize; iSite++) {
            for (int iMat = 0; iMat < 18; iMat++) {
                if (lattice[mu][iSite][iMat] == 0 || lattice[mu][iSite][iMat] == 1) {
                    printf("\n\nERROR!\n\n"); // Sanity check for us loading correct matrix...
                }
            }
        }
    }

    ///////////////////////////
    //////// PLAQUETTE ////////
    ///////////////////////////

    // Initializes samples for the
    Lattice <SU3> Temp1(dim), Temp2(dim);
    Temp1.zeros();
    Temp2.zeros();

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu+1; nu < 4; nu++) {
            Temp1 = lattice[mu];
            Temp1 *= shift(lattice[nu],FORWARDS,mu);
            Temp1 *= shift(lattice[mu],FORWARDS,nu).inv();
            Temp1 *= lattice[nu].inv();
            Temp2 += Temp1;
        }
    }
    double observable = sumRealTrace(Temp2)/m_multiplicationFactor;
    MPI_Allreduce(&observable,&observable,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    observable /= double(Parallel::Communicator::getNumProc());
    printf("\nPlaquette = %20.16f\n",observable);

    ///////////////////////////
    //// TOPOLOGICAL CHARGE ///
    ///////////////////////////
    double tempDiag1, tempDiag2;
    int mu = 0;
    Lattice<SU3> clov1(dim), clov2(dim), U2Temp(dim), U3Temp(dim);
    for (int nu = 1; nu < 4; nu++)
    {
        // First clover. Definition from R Wohler 1985, more symmetric than other methods.
        // First leaf
        Temp1 = lattice[mu];
        Temp1 *= shift(lattice[nu],FORWARDS,mu);
        Temp1 *= shift(lattice[mu],FORWARDS,nu).inv();
        Temp1 *= lattice[nu].inv();
        Temp2 = Temp1;
//        U1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
//        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
//        U1 *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
//        U1 *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
//        clov1 = U1;

        // Adds plaquette leaf
//        m_plaquettes[2*(nu-1)] = U1;

        // Retrieves beforehand in order to reduce number of communications by 2.
        U2Temp = shift(lattice[nu],BACKWARDS,nu);
        U3Temp = shift(lattice[mu],BACKWARDS,mu).inv();

        // Second leaf
        Temp1 = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
        Temp1 *= Parallel::Communicator::getNeighboursNeighbourLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
        Temp1 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,mu).inv();
//        U2 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
        Temp1 *= U2Temp;
        Temp2 -= Temp1;


        // Third leaf
//        U3 = Parallel::Communicator::getNegativeLink(lattice,m_position,mu,muIndex,mu).inv();
        Temp1 = U3Temp;
        Temp1 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,nu).inv();
        Temp1 *= Parallel::Communicator::getNeighboursNeighbourNegativeLink(lattice,m_position,mu,muIndex,nu,nuIndex,mu);
//        U3 *= Parallel::Communicator::getNegativeLink(lattice,m_position,nu,nuIndex,nu);
        Temp1 *= U2Temp;
        Temp2 += U3;

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

    ///////////////////////////
    ///////// ENERGY //////////
    ///////////////////////////

}
