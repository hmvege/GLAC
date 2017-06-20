#include "action.h"
#include "links.h"
#include "su3.h"
#include "functions.h"

// TEMP
#include <iostream>
using std::cout;
using std::endl;

//Action::Action()
//{
////    muIndex = new int[4];
////    nuIndex = new int[4];
////    for (int i = 0; i < 4; i++) {
////        muIndex[i] = 0;
////        nuIndex[i] = 0;
////    }
//}

Action::~Action()
{
//    delete [] muIndex;
//    delete [] nuIndex;
}

Action::Action(int N)
{
    m_N = N;
//    beta = new_beta;
//    muIndex = new int[4];
//    nuIndex = new int[4];
//    for (int i = 0; i < 4; i++) {
//        muIndex[i] = 0;
//        nuIndex[i] = 0;
//    }
}

double Action::getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu)
{
    cout << "If you are seeing this, something is wrong!" << endl;
    exit(1);
    /*
     * Takes the entire lattice and the updated spacetime index
     * Arguments:
     *  lattice : Lattice of N^4 * 4 SU3 matrices
     *  U       : Updated SU3 matrix
     *  n       : spacetime index for updated SU3 matrix
     *  mu      : Lorentz index for updated SU3 matrix
     */
//    double S = 0;
//    SU3 A;
//    SU3 tr;
//    for (int nu = 0; nu < 4; nu++)
//    {
//        if (mu == nu) continue;
//        lorentzIndex(mu,muIndex);
//        lorentzIndex(nu,nuIndex);

//        A += lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], N)].U[nu]
//                *inverse(lattice[stapleIndex(i+muIndex[0]+nuIndex[0],j+muIndex[1]+nuIndex[1],k+muIndex[2]+nuIndex[2],l+muIndex[3]+nuIndex[3], N)].U[mu])
//                *inverse(lattice[stapleIndex(i+nuIndex[0],j+nuIndex[1],k+nuIndex[2],l+nuIndex[3], N)].U[nu])
//                + inverse(lattice[stapleIndex(i+muIndex[0],j+muIndex[1],k+muIndex[2],l+muIndex[3], N)].U[nu])
//                *inverse(lattice[stapleIndex(i+muIndex[0]-nuIndex[0],j+muIndex[1]-nuIndex[1],k+muIndex[2]-nuIndex[2],l+muIndex[3]-nuIndex[3], N)].U[mu])
//                *lattice[stapleIndex(i-nuIndex[0],j-nuIndex[1],k-nuIndex[2],l-nuIndex[3], N)].U[nu];
//    }
//    tr = (U - lattice[index(i,j,k,l,N)].U[mu])*A;
//    for (int i = 0; i < 3; i++)
//    {
//        S += tr.mat[i*3+i].re;
//    }
//    return beta/3.0*S; // Should be N=3 as in the Gauge symmetry?
    return 1.0;
}
