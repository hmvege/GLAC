#include "mastersampler.h"

#include "parallelization/communicator.h"

#include "io/fieldio.h"

//using namespace LatticeOperations;

MasterSampler::MasterSampler()
{
    m_multiplicationFactor = 1;
//    int Nx = 2, Ny = 1, Nz = 1, Nt = 1;
//    std::vector<int> dim = {Nx,Ny,Nz,Nt};
//    int latticeSize = Nx*Ny*Nz*Nt;
//    Lattice <SU3> A, B;
//    A.allocate(dim);
//    B.allocate(dim);
//    A.identity();
//    B.identity();
//    printf("\nA.N=%d B.N=%d\n",A.m_latticeSize,B.m_latticeSize);
//    B *= 1000;
//    B[latticeSize-1][1] = 10;
//    B[latticeSize-1][17] = 10;
//    Parallel::Communicator::setBarrier();
//    if (Parallel::Communicator::getProcessRank() == 0) A[latticeSize-1].print();
//    A += B;
//    Parallel::Communicator::setBarrier();
//    if (Parallel::Communicator::getProcessRank() == 0) A[latticeSize-1].print();
//    A -= B;
//    Parallel::Communicator::setBarrier();
//    if (Parallel::Communicator::getProcessRank() == 0) A[latticeSize-1].print();
//    A = A + B;

//    for (int i = 0; i < latticeSize; i++) {
//        A[i] *= 2;
//    }
//    printf("\nSumming!");
//    SU3 a = sum(A);
//    a.print();
////    A[0].print();
//    printf("\nShifting!");
//    Parallel::Communicator::setBarrier();
//    printf("\nShifting!");
//    Parallel::Communicator::setBarrier();
//    int mu = 0;
//    B.zeros();
//    A.identity();
//    A = shift(B,FORWARDS,mu);
//    B[0].print();
//    B[1].print();
//    A[0].print();
//    A[1].print();

}

void MasterSampler::calculate()
{
    Lattice <SU3> lattice[4];
    IO::FieldIO::loadLatticeFieldConfiguration("",lattice);
    int Nx = 2, Ny = 1, Nz = 1, Nt = 1;
    std::vector<int> dim = {Nx,Ny,Nz,Nt};
    Lattice <SU3> PTemp(dim), P(dim);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu+1; nu < 4; nu++) {
            PTemp = lattice[mu];
            PTemp *= shift(lattice[nu],FORWARDS,mu);
            PTemp *= shift(lattice[mu],FORWARDS,nu).inv();
            PTemp *= lattice[nu].inv();
            P += PTemp;
        }
    }
    double observable = sum(realTrace(P))/m_multiplicationFactor;
    printf("\nPlaquette = %20.16f\n",observable);
}
