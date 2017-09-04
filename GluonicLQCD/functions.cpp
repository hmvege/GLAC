#include "matrices/su3.h"
#include "matrices/su2.h"

//int index(int i, int j, int k, int l, int N, int N_T)
//{
//    /*
//     * Function for contigious memory allocation.
//     */
////    if ((i < 0) || (j < 0) || (k < 0) || (l < 0) || (i > N) || (j > N) || (k > N) || (l > N)) { //TEST FOR ENSURING CORRECT INDICES
////        std::cout << "Index error!" << std::endl;
////        exit(1);
////    }
//    int returnVal = (N_T*(N*(N*i + j) + k) + l);

////    if (returnVal == (N*N*N*2*N)) {
////        std::cout << "good" << std::endl;
////        exit(1);
////    }
//    return returnVal;
//}

int getIndex(int i, int j, int k, int l, int Ny, int Nz, int Nt)
{
    /*
     * N[0] = x-direction
     * N[1] = y-direction
     * N[2] = z-direction
     * N[3] = t-direction
     * Function for contigious memory allocation for uneven dimensions
     */
    return (Nt*(Nz*(Ny*i + j) + k) + l);
}

SU3 inverse(SU3 U)
{
    /*
     * Gets the inverse of a SU3 matrix, MAKE CLASS BASED?!
     */
    SU3 UInverse;
    UInverse.copy(U);
    UInverse.conjugate();
    UInverse.transpose();
    return UInverse;
}

SU2 SU2Inverse(SU2 U)
{
    SU2 UInverse;
    UInverse.copy(U);
    UInverse.conjugate();
    UInverse.transpose();
    return UInverse;
}

int stapleIndex(int i, int j, int k, int l, int *N)
{
    /*
     * Unit vector for lorentz indexes.
     * N[0] = x-direction
     * N[1] = y-direction
     * N[2] = z-direction
     * N[3] = t-direction
     */
//    return getIndex((i+N[0]-2) % (N[0]-2), (j+N[1]-2) % (N[1]-2), (k+N[2]-2) % (N[2]-2), (l+N[3]-2) % (N[3]-2), N[1], N[2], N[3]);
    return getIndex((i+N[0]) % N[0], (j+N[1]) % N[1], (k+N[2]) % N[2], (l+N[3]) % N[3], N[1], N[2], N[3]);
}

void lorentzIndex(int mu, int *lorentzIndices)
{
    /*
     * Fills the lorentz array with correct indices
     */
    for (int i = 0; i < 4; i++)
    {
        if (i==mu)
        {
            lorentzIndices[i] = 1;
        }
        else
        {
            lorentzIndices[i] = 0;
        }
    }
}

complex SU3Determinant(SU3 U)
{
    complex det;
    det += U[2]*(U[3]*U[7] - U[6]*U[4]);
    det += U[5]*(U[6]*U[1] - U[0]*U[7]);
    det += U[8]*(U[0]*U[4] - U[3]*U[1]);
    return det;
}

bool compareSU3(SU3 A, SU3 B)
{
    for (int i = 0; i < 9; i++) {
        if ((A.mat[i].re() != B.mat[i].re()) || (A.mat[i].im() != B.mat[i].im())) {
            return false;
        }
    }
    return true;
}
