#include "matrices/su3.h"
#include "matrices/su2.h"

int index(int i, int j, int k, int l, int N, int N_T)
{
    /*
     * Function for contigious memory allocation.
     */
//    if ((i < 0) || (j < 0) || (k < 0) || (l < 0) || (i > N) || (j > N) || (k > N) || (l > N)) { //TEST FOR ENSURING CORRECT INDICES
//        std::cout << "Index error!" << std::endl;
//        exit(1);
//    }
    int returnVal = (N_T*(N*(N*i + j) + k) + l);

//    if (returnVal == (N*N*N*2*N)) {
//        std::cout << "good" << std::endl;
//        exit(1);
//    }
    return returnVal;
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

int stapleIndex(int i, int j, int k, int l, int N, int N_T)
{
    /*
     * Unit vector for lorentz indexes
     */
    return index((i+N) % N, (j+N) % N, (k+N) % N, (l+N_T) % N_T, N, N_T);;
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
