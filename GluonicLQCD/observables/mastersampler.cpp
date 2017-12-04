#include "mastersampler.h"

#include <cmath>
#include "parallelization/communicator.h"
#include "config/parameters.h"
#include "io/fieldio.h"

//using namespace LatticeOperations;

MasterSampler::MasterSampler()
{
    m_latticeSize = Parameters::getLatticeSize();
    m_plaqMultiplicationFactor = 1.0/(18.0*double(m_latticeSize));
    m_topcMultiplicationFactor = 1.0/(16*16*M_PI*M_PI);
    m_energyMultiplicationFactor = 1.0/double(m_latticeSize);

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
//    m_multiplicationFactor = 18.0*double(Parameters::getSubLatticeSize());
    IO::FieldIO::init();
    std::vector<int> dim = {int(N[0]),int(N[1]),int(N[2]),int(N[3])};
    printf("\n%d %d %d %d \n",int(N[0]),int(N[1]),int(N[2]),int(N[3]));
    for (int mu = 0; mu < 4; mu++) {
        lattice[mu].allocate(dim);
        lattice[mu].identity();
    }

    m_latticeSize = Parameters::getSubLatticeSize();
    m_plaqMultiplicationFactor = 1.0/(18.0*double(m_latticeSize));
    m_topcMultiplicationFactor = 1.0/(16*16*M_PI*M_PI);
    m_energyMultiplicationFactor = 1.0/double(m_latticeSize);

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
    //// SYMMETRIC CLOVER /////
    ///////////////////////////
    double topCharge = 0;
    double energy = 0;
    double plaquette = 0;
    Lattice <double> tempDiag1;
    int mu = 0;
    Lattice<SU3> clov1(dim), clov2(dim), U2Temp(dim), U3Temp(dim), Temp1(dim);

    for (int nu = 1; nu < 4; nu++)
    {
        // First clover. Definition from R Wohler 1985, more symmetric than other methods.
        // First leaf
        Temp1 = lattice[mu];
        Temp1 *= shift(lattice[nu],FORWARDS,mu);
        Temp1 *= shift(lattice[mu],FORWARDS,nu).inv();
        Temp1 *= lattice[nu].inv();
        clov1 = Temp1;

        // Adds plaquette leaf
        plaquette += sumRealTrace(clov1);

        // Retrieves beforehand in order to reduce number of communications by 2.
        U2Temp = shift(lattice[nu],BACKWARDS,nu);
        U3Temp = shift(lattice[mu],BACKWARDS,mu).inv();

        // Second leaf
        Temp1 = lattice[mu];
        Temp1 *= shift(shift(lattice[nu],FORWARDS,mu),BACKWARDS,nu).inv();
        Temp1 *= shift(lattice[mu],BACKWARDS,nu).inv();
        Temp1 *= U2Temp;
        clov1 -= Temp1;

        // Third leaf
        Temp1 = U3Temp;
        Temp1 *= shift(shift(lattice[nu],BACKWARDS,mu),BACKWARDS,mu).inv();
        Temp1 *= shift(shift(lattice[mu],BACKWARDS,mu),BACKWARDS,nu);
        Temp1 *= U2Temp;
        clov1 += Temp1;

        // Fourth leaf
        Temp1 = U3Temp;
        Temp1 *= shift(lattice[nu],BACKWARDS,mu);
        Temp1 *= shift(shift(lattice[mu],FORWARDS,nu),BACKWARDS,mu);
        Temp1 *= lattice[nu].inv();
        clov1 -= Temp1;

        int rho = nu % 3;
        rho++;
        int sigma = rho % 3;
        sigma++;

        // Second clover
        // First leaf
        Temp1 = lattice[rho];
        Temp1 *= shift(lattice[sigma],FORWARDS,rho);
        Temp1 *= shift(lattice[rho],FORWARDS,sigma).inv();
        Temp1 *= lattice[sigma].inv();
        clov2 = Temp1;

        // Gets lattice for temp use
        U2Temp = shift(lattice[sigma],BACKWARDS,sigma);
        U3Temp = shift(lattice[rho],BACKWARDS,rho).inv();

        // Adds another leaf to the plaquette
        plaquette += sumRealTrace(clov2);

        // Second leaf
        Temp1 = lattice[rho];
        Temp1 *= shift(shift(lattice[sigma],FORWARDS,rho),BACKWARDS,sigma).inv();
        Temp1 *= shift(lattice[rho],BACKWARDS,sigma).inv();
        Temp1 *= U2Temp;
        clov2 -= Temp1;

        // Third leaf
        Temp1 = U3Temp;
        Temp1 *= shift(shift(lattice[sigma],BACKWARDS,rho),BACKWARDS,sigma).inv();
        Temp1 *= shift(shift(lattice[rho],BACKWARDS,rho),BACKWARDS,sigma);
        Temp1 *= U2Temp;
        clov2 += Temp1;

        // Fourth leaf
        Temp1 = U3Temp;
        Temp1 *= shift(lattice[sigma],BACKWARDS,rho);
        Temp1 *= shift(shift(lattice[rho],FORWARDS,sigma),BACKWARDS,rho);
        Temp1 *= lattice[sigma].inv();
        clov2 += Temp1;

        // Makes first clover anti hermitian and traceless
        Temp1 = clov1.inv();
        clov1 -= Temp1;
        tempDiag1 = imagTrace(clov1)/3.0;
        clov1 = subtractImag(clov1,tempDiag1);

        // Makes second clover anti hermitian and traceless
        Temp1 = clov2.inv();
        clov2 -= Temp1;
        tempDiag1 = imagTrace(clov2)/3.0;
        clov2 = subtractImag(clov2,tempDiag1);

        // Picks up the topological charge
        topCharge -= sumRealTraceMultiplication(clov1,clov2);

        // Picks up the action density
        energy += sumRealTraceMultiplication(clov1,clov1);
        energy += sumRealTraceMultiplication(clov2,clov2);
    }

    ///////////////////////////
    //////// PLAQUETTE ////////
    ///////////////////////////

//    // Initializes samples for the
//    Lattice<SU3>Temp2(dim);
//    Temp1.zeros();
//    Temp2.zeros();

//    for (int mu = 0; mu < 4; mu++) {
//        for (int nu = mu+1; nu < 4; nu++) {
//            Temp1 = lattice[mu];
//            Temp1 *= shift(lattice[nu],FORWARDS,mu);
//            Temp1 *= shift(lattice[mu],FORWARDS,nu).inv();
//            Temp1 *= lattice[nu].inv();
//            Temp2 += Temp1;
//        }
//    }
//    double observable = sumRealTrace(Temp2)*m_plaqMultiplicationFactor;

//    MPI_Allreduce(&observable,&observable,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    observable /= double(Parallel::Communicator::getNumProc());
//    if (Parallel::Communicator::getProcessRank() == 0) printf("\nPlaquette = %20.16f\n",observable);

//    plaquette *= m_plaqMultiplicationFactor;
//    MPI_Allreduce(&observable,&observable,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    observable /= double(Parallel::Communicator::getNumProc());
//    printf("\nPlaquette = %20.16f\n",observable);
    plaquette *= m_plaqMultiplicationFactor;
    MPI_Allreduce(&plaquette,&plaquette,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    plaquette /= double(Parallel::Communicator::getNumProc());
    if (Parallel::Communicator::getProcessRank() == 0) printf("\nPlaquette           = %20.16f",plaquette);

    ///////////////////////////
    //// TOPOLOGICAL CHARGE ///
    ///////////////////////////
    topCharge *= m_topcMultiplicationFactor;
    MPI_Allreduce(&topCharge,&topCharge,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    topCharge /= double(Parallel::Communicator::getNumProc());
    if (Parallel::Communicator::getProcessRank() == 0) printf("\nTopological charge  = %20.16f",topCharge);

    ///////////////////////////
    ///////// ENERGY //////////
    ///////////////////////////
    energy *= m_energyMultiplicationFactor;
    MPI_Allreduce(&energy,&energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    energy /= double(Parallel::Communicator::getNumProc());
    if (Parallel::Communicator::getProcessRank() == 0) printf("\nEnergy              = %20.16f",energy);

}
