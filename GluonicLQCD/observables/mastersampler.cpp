#include "mastersampler.h"


MasterSampler::MasterSampler()
{
    m_multiplicationFactor = 1;
    int Nx = 2, Ny = 1, Nz = 1, Nt = 1;
    int N = Nx*Ny*Nz*Nt;
    Lattice <SU3> A, B;
    A.allocate(N,Nx,Ny,Nz,Nt);
    B.allocate(N,Nx,Ny,Nz,Nt);
    for (int i = 0; i < N; i++) {
        A[i].identity();
        B[i].identity();
    }
    printf("\nA.N=%d B.N=%d",A.m_N,B.m_N);
    B *= 1000;
    B[N-1][1] = 10;
    B[N-1][17] = 10;
    A += B;
//    A[0].print();
    A[N-1].print();
//    for (int i = 0; i < N; i++) {
//        A[i].print();
//    }
//    exit(1);
    A = A + B;
    exit(1);
    for (int i = 0; i < N; i++) {
        A[i] *= 2;
    }
    A[0].print();
    B = A;
    B[0].print();

}

void MasterSampler::calculate()
{
//    P.zeros();
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
