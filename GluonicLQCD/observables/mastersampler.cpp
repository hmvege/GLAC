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

    // Initializes samples for the
    Lattice <SU3> PTemp(dim), P(dim);
    P.zeros();
    PTemp.zeros();

    P.zeros();
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu+1; nu < 4; nu++) {
            PTemp = lattice[mu];
            PTemp *= shift(lattice[nu],FORWARDS,mu);
            PTemp *= shift(lattice[mu],FORWARDS,nu).inv();
            PTemp *= lattice[nu].inv();
            P += PTemp;
        }
    }
    double observable = sumRealTrace(P)/m_multiplicationFactor;
    MPI_Allreduce(&observable,&observable,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    observable /= double(Parallel::Communicator::getNumProc());
    printf("\nPlaquette = %20.16f\n",observable);
}
