#include "mastersampler.h"


MasterSampler::MasterSampler()
{
    m_multiplicationFactor = 1;
    int Nx = 2, Ny = 1, Nz = 1, Nt = 1;
    std::vector<int> dim = {Nx,Ny,Nz,Nt};
    int latticeSize = Nx*Ny*Nz*Nt;
    Lattice <SU3> A, B;
    A.allocate(dim);
    B.allocate(dim);
    for (int i = 0; i < latticeSize; i++) {
        A[i].identity();
        B[i].identity();
    }
    printf("\nA.N=%d B.N=%d",A.m_latticeSize,B.m_latticeSize);
    B *= 1000;
    B[latticeSize-1][1] = 10;
    B[latticeSize-1][17] = 10;
    A += B;
    A += B;
//    for (int i = 0; i < N; i++) {
//        A[i].print();
//    }
    A = A + B;
//    A[0].print();
//    A[N-1].print();
    for (int i = 0; i < latticeSize; i++) {
        A[i] *= 2;
    }
    A.sum();
    A[0].print();
    B = A;
    B[0].print();

}

void MasterSampler::calculate()
{
    int Nx = 2, Ny = 1, Nz = 1, Nt = 1;
    std::vector<int> dim = {Nx,Ny,Nz,Nt};
//    int N = Nx*Ny*Nz*Nt;
    Lattice <SU3> A(dim);
//    for (int mu = 0; mu < 4; mu++) {
//        updateMuIndex(mu); // Inline function
//        for (int nu = mu+1; nu < 4; nu++) {
//            updateNuIndex(nu); // Inline function
//            PTemp = lattice[Parallel::Index::getIndex(i,j,k,l)].U[mu];
//            PTemp *= Parallel::Communicator::getPositiveLink(lattice,m_position,mu,muIndex,nu);
//            PTemp *= Parallel::Communicator::getPositiveLink(lattice,m_position,nu,nuIndex,mu).inv();
//            PTemp *= lattice[Parallel::Index::getIndex(i,j,k,l)].U[nu].inv();
//            P += PTemp;
//        }
//    }
//    m_observable->m_observables[iObs] = (P.mat[0] + P.mat[8] + P.mat[16])/m_multiplicationFactor;
}
