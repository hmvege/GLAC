#ifndef TESTCORE_H
#define TESTCORE_H

#include "math/lattice.h"
#include "math/latticemath.h"
#include "math/matrices/su3matrixgenerator.h"

class TestCore
{
protected:
    // SU3 Test variables
    SU3 U1, U2, U3, UAdd, USub, UMul, UC, UT, UCT, UTrace, U_RST, U1Hermitian, U1AntiHermitian;

    // SU2 Test variables
    SU2 s1, s2, s3, sAdd, sSub, sMul, sC, sT, sCT, s_r, s_s, s_t;

    // Complex test variables
    complex z1, z2, zAdd, zSub, zMul, zDiv, zConj, zSetToMinus;
    double zNorm, zNormSquared;

    // Lattice test variables. Since only the matrix changes, we can reuse previous variables
    std::vector<unsigned int> m_dim;
    Lattice<SU3> latticeSU3_U1, latticeSU3_U2;
    Lattice<complex> latticeComplex_z1, latticeComplex_z2;
    double latticeDoubleValue1, latticeDoubleValue2;
    Lattice<double> latticeDouble1, latticeDouble2;

    // Limit we demand matrix properties to be retained at
    double m_eps = 2*1e-14;
    double m_tracedMatrix;

    // Machine precision
    double m_machineEpsilon = 1e-16;

    // Variables used for full lattice testing
    int m_processRank;
    int m_numprocs;
    unsigned int m_N, m_NT;
    unsigned long int m_subLatticeSize;

    // Verbose storage
    bool m_verbose = false;

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;

    // Basic matrix property testers
    bool operationSU2Test(SU2 results, SU2 solution, std::string operation);
    bool operationSU3Test(SU3 results, SU3 solution, std::string operation);

    // Inline matrix comparing functions
    inline bool compareSU3(const SU3 A, const SU3 B);
    inline bool compareSU2(const SU2 A, const SU2 B);
    inline bool compareComplex(const complex a, const complex b);

    // Inline complex dot product function
    inline complex dot(complex * a, complex * b);
    inline complex dot2(complex * a, complex * b);

public:
    TestCore();
    ~TestCore();
};

////////////////////////////////////
////// Comparison functions ////////
////////////////////////////////////
inline bool TestCore::compareComplex(const complex a, const complex b)
{
    /*
     * Function that compares two complex numbers and returns false if they are different.
     */
    for (int i = 0; i < 2; i++)
    {
        if (fabs(a.z[i] - b.z[i]) > m_machineEpsilon)
        {
            return false;
        }
    }
    return true;
}

inline bool TestCore::compareSU2(const SU2 A, const SU2 B)
{
    /*
     * Function that compares two SU2 matrices and returns false if they are not EXACTLY the same.
     */
    for (int i = 0; i < 8; i++)
    {
        if (fabs(A.mat[i] - B.mat[i]) > m_machineEpsilon)
        {
            return false;
        }
    }
    return true;
}

inline bool TestCore::compareSU3(const SU3 A, const SU3 B)
{
    /*
     * Function that compares two SU3 matrices and returns false if they are not EXACTLY the same.
     */
    for (int i = 0; i < 18; i++)
    {
        if (fabs(A.mat[i] - B.mat[i]) > m_machineEpsilon)
        {
            return false;
        }
    }
    return true;
}

////////////////////////////////////
/////// Vector dot prodcts /////////
////////////////////////////////////
// Inline dot product between two columns of a 3x3 matrix
inline complex TestCore::dot(complex * a, complex * b)
{
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex returnSum(0,0);
    for (int i = 0; i < 3; i++)
    {
        returnSum += a[i].c()*b[i];
    }
    return returnSum;
}

// Inline dot product between two columns of a 2x2 matrix
inline complex TestCore::dot2(complex * a, complex * b)
{
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex returnSum(0,0);
    for (int i = 0; i < 2; i++)
    {
        returnSum += a[i].c()*b[i];
    }
    return returnSum;
}

#endif // TESTCORE_H
