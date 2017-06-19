#include "action.h"
#include "links.h"
#include "su3.h"
#include "functions.h"

Action::Action()
{

}

Action::Action(int latticeSize, double new_a, double new_beta)
{
    N = latticeSize;
    a = new_a;
    beta = new_beta;
}

double Action::getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu)
{
    /*
     * Takes the entire lattice and the updated spacetime index
     * Arguments:
     *  lattice : Lattice of N^4 * 4 SU3 matrices
     *  U       : Updated SU3 matrix
     *  n       : spacetime index for updated SU3 matrix
     *  mu      : Lorentz index for updated SU3 matrix
     */
    double S;
    SU3 A;
    SU3 tr;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu) continue;
//        A += lattice[index(i,j,k,l,N)].U[nu]*inverse(lattice[index(i,j,k,l, N)].U[mu])*inverse(lattice[index(i,j,k,l, N)].U[nu])
//                + inverse(lattice[index(i,j,k,l, N)].U[nu])*inverse(lattice[index(i,j,k,l, N)].U[mu])*lattice[index(i,j,k,l, N)].U[nu];
        A += lattice[index(i,j,k,l,N)].U[nu];
        FORTSETTE HER: hvordan mappe mu,nu positive og negative til korrekte indexer med periodiske boundary conditions?
    }
    tr = (U - lattice[index(i,j,k,l,N)].U[mu])*A;
    for (int i = 0; i < 3; i++)
    {
        S += tr.mat[i*3+i];
    }
    return beta/double(N)*S; // Correct N?!?!?!? Think it should be 6, the number of staples
}

int Action::stapleIndex(int i, int j, int k, int l, int N, int mu)
{
    /*
     * Unit vector for lorentz indexes
     */
    if (mu==0) {
        return index((i+1) % N,j,k,l,N);
    } else if (mu==1) {
        return index(i,(j+1) % N,k,l,N);
    } else if (mu==2) {
        return index(i,j,(k+1) % N,l,N);
    } else {
        return index(i,j,k,(l+1) % N,N);
    }
}
