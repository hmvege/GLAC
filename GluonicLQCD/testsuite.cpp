#include "testsuite.h"

using std::cout;
using std::endl;

TestSuite::TestSuite(SU3MatrixGenerator *SU3Gen)
{
    /*
     * Class for running unit tests.
     */
    // Initiating the Mersenne Twister random number generator
    std::mt19937_64 gen(std::time(nullptr));
    std::uniform_real_distribution<double> uni_dist(0,1);
    m_generator = gen;
    m_uniform_distribution = uni_dist;

    // Initiating the SU3 Matrix generator
//    SU3MatrixGenerator SU3Gen(0.24, std::time(nullptr)+1);
    m_SU3Generator = SU3Gen;

    //// 3x3 COMPLEX MATRICES
    // Setting up matrix U1
    U1.setComplex(complex(1,1),0);
    U1.setComplex(complex(1,2),2);
    U1.setComplex(complex(1,3),4);
    U1.setComplex(complex(2,1),6);
    U1.setComplex(complex(2,2),8);
    U1.setComplex(complex(2,3),10);
    U1.setComplex(complex(3,1),12);
    U1.setComplex(complex(3,2),14);
    U1.setComplex(complex(3,3),16);
    // Setting up matrix U2
    U2.setComplex(complex(4,4),0);
    U2.setComplex(complex(4,5),2);
    U2.setComplex(complex(4,6),4);
    U2.setComplex(complex(5,4),6);
    U2.setComplex(complex(5,5),8);
    U2.setComplex(complex(5,6),10);
    U2.setComplex(complex(6,4),12);
    U2.setComplex(complex(6,5),14);
    U2.setComplex(complex(6,6),16);
    // Addition
    UAdd.setComplex(complex(5,5),0);
    UAdd.setComplex(complex(5,7),2);
    UAdd.setComplex(complex(5,9),4);
    UAdd.setComplex(complex(7,5),6);
    UAdd.setComplex(complex(7,7),8);
    UAdd.setComplex(complex(7,9),10);
    UAdd.setComplex(complex(9,5),12);
    UAdd.setComplex(complex(9,7),14);
    UAdd.setComplex(complex(9,9),16);
    // Subtraction
    USub.setComplex(complex(-3,-3),0);
    USub.setComplex(complex(-3,-3),6);
    USub.setComplex(complex(-3,-3),12);
    USub.setComplex(complex(-3,-3),2);
    USub.setComplex(complex(-3,-3),8);
    USub.setComplex(complex(-3,-3),14);
    USub.setComplex(complex(-3,-3),4);
    USub.setComplex(complex(-3,-3),10);
    USub.setComplex(complex(-3,-3),16);
    // Multiplication
    UMul.setComplex(complex(-9,44),0);
    UMul.setComplex(complex(-15,47),2);
    UMul.setComplex(complex(-21,50),4);
    UMul.setComplex(complex(6,56),6);
    UMul.setComplex(complex(0,62),8);
    UMul.setComplex(complex(-6,68),10);
    UMul.setComplex(complex(21,68),12);
    UMul.setComplex(complex(15,77),14);
    UMul.setComplex(complex(9,86),16);
    // Conjugation
    UC.setComplex(complex(1,-1),0);
    UC.setComplex(complex(1,-2),2);
    UC.setComplex(complex(1,-3),4);
    UC.setComplex(complex(2,-1),6);
    UC.setComplex(complex(2,-2),8);
    UC.setComplex(complex(2,-3),10);
    UC.setComplex(complex(3,-1),12);
    UC.setComplex(complex(3,-2),14);
    UC.setComplex(complex(3,-3),16);
    // Transpose
    UT.setComplex(complex(1,1),0);
    UT.setComplex(complex(2,1),2);
    UT.setComplex(complex(3,1),4);
    UT.setComplex(complex(1,2),6);
    UT.setComplex(complex(2,2),8);
    UT.setComplex(complex(3,2),10);
    UT.setComplex(complex(1,3),12);
    UT.setComplex(complex(2,3),14);
    UT.setComplex(complex(3,3),16);
    // Conjugate transpose
    UCT.setComplex(complex(1,-1),0);
    UCT.setComplex(complex(2,-1),2);
    UCT.setComplex(complex(3,-1),4);
    UCT.setComplex(complex(1,-2),6);
    UCT.setComplex(complex(2,-2),8);
    UCT.setComplex(complex(3,-2),10);
    UCT.setComplex(complex(1,-3),12);
    UCT.setComplex(complex(2,-3),14);
    UCT.setComplex(complex(3,-3),16);
    /// 2x2 COMPLEX MATRICE
    // s1
    s1.setComplex(complex(1,1),0);
    s1.setComplex(complex(1,2),2);
    s1.setComplex(complex(2,1),4);
    s1.setComplex(complex(2,2),6);
    // s2
    s2.setComplex(complex(4,4),0);
    s2.setComplex(complex(4,5),2);
    s2.setComplex(complex(5,4),4);
    s2.setComplex(complex(5,5),6);
    // Adding
    sAdd.setComplex(complex(5,5),0);
    sAdd.setComplex(complex(5,7),2);
    sAdd.setComplex(complex(7,5),4);
    sAdd.setComplex(complex(7,7),6);
    // Subtracting
    sSub.setComplex(complex(-3,-3),0);
    sSub.setComplex(complex(-3,-3),2);
    sSub.setComplex(complex(-3,-3),4);
    sSub.setComplex(complex(-3,-3),6);
    // Multiplying
    sMul.setComplex(complex(-3,22),0);
    sMul.setComplex(complex(-6,24),2);
    sMul.setComplex(complex(6,30),4);
    sMul.setComplex(complex(3,34),6);
    // Conjugate of s1
    sC.setComplex(complex(1,-1),0);
    sC.setComplex(complex(1,-2),2);
    sC.setComplex(complex(2,-1),4);
    sC.setComplex(complex(2,-2),6);
    // Transpose of s1
    sT.setComplex(complex(1,1),0);
    sT.setComplex(complex(2,1),2);
    sT.setComplex(complex(1,2),4);
    sT.setComplex(complex(2,2),6);
    // Complex conjugate of s1
    sCT.setComplex(complex(1,-1),0);
    sCT.setComplex(complex(2,-1),2);
    sCT.setComplex(complex(1,-2),4);
    sCT.setComplex(complex(2,-2),6);
}

void TestSuite::runFullTestSuite(bool verbose)
{
    if (verbose)
    {
        // Original matrices
        cout << "Original matrices\nU1 =" << endl;
        U1.print();
        cout << "U2 =" << endl;
        U2.print();
    }

    if (run3x3MatrixTests(verbose) & run2x2MatrixTests(verbose) & runSU2Tests(verbose) & runSU3Tests(verbose)) {
        cout << "SUCCESS: All tests passed." << endl;
    } else {
        cout << "FAILURE: One or more test(s) failed." << endl;
    }
}

bool TestSuite::runSU2Tests(bool verbose)
{
    cout << "Running SU2 property tests." << endl;
    bool passed = (testSU2Hermicity(verbose) && testSU2Orthogonality(verbose) && testSU2Norm(verbose)
                   && testSU2Determinant(verbose));
    if (passed) {
        cout << "PASSED: SU2 properties." << endl;
    } else {
        cout << "FAILED: SU2 properties." << endl;
    }

    return passed;
}

bool TestSuite::runSU3Tests(bool verbose)
{
    cout << "Running SU3 property tests." << endl;
    bool passed = (testSU3Hermicity(verbose) && testSU3Orthogonality(verbose) && testSU3Norm(verbose)
                   && testSU3TraceMultiplication(verbose) && testSU3Determinant(verbose));
    if (passed) {
        cout << "PASSED: SU3 properties." << endl;
    } else {
        cout << "FAILED: SU3 properties." << endl;
    }

    return passed;
}

bool TestSuite::run2x2MatrixTests(bool verbose)
{
    cout << "Running basic SU2 matrix tests." << endl;
    bool passed = (testSU2Addition(verbose) && testSU2Subtraction(verbose) && testSU2Multiplication(verbose)
                   && testSU2Conjugation(verbose) && testSU2Transpose(verbose));
    if (passed) {
        cout << "PASSED: Basic SU2 matrix operations." << endl;
    } else {
        cout << "FAILED: Basic SU2 matrix operations." << endl;
    }
    return passed;
}

bool TestSuite::run3x3MatrixTests(bool verbose)
{
    cout << "Running basic SU3 matrix tests." << endl;
    bool passed = (testSU3Addition(verbose) && testSU3Subtraction(verbose) && testSU3Multiplication(verbose)
                   && testSU3Conjugation(verbose) && testSU3Transpose(verbose));
    if (passed) {
        cout << "PASSED: Basic SU3 matrix operations." << endl;
    } else {
        cout << "FAILED: Basic SU3 matrix operations." << endl;
    }
    return passed;
}

// Comparison functions ===================================
inline bool TestSuite::compareSU2(SU2 A, SU2 B)
{
    /*
     * Function that compares two SU2 matrices and returns false if they are not EXACTLY the same.
     */
    for (int i = 0; i < 8; i++) {
        if (A.mat[i] != B.mat[i]) {
            return false;
        }
    }
    return true;
}

inline bool TestSuite::compareSU3(SU3 A, SU3 B)
{
    /*
     * Function that compares two SU3 matrices and returns false if they are not EXACTLY the same.
     */
    for (int i = 0; i < 18; i++) {
        if (A.mat[i] != B.mat[i]) {
            return false;
        }
    }
    return true;
}

// Inline dot product between two columns
inline complex dot(complex * a, complex * b) {
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex returnSum(0,0);
    for (int i = 0; i < 3; i++) {
        returnSum += a[i]*b[i].c();
    }
    return returnSum;
}

// Matrix operation tester ================================
bool TestSuite::operationSU2Test(bool verbose, SU2 results, SU2 solution, std::string operation)
{
    if (compareSU2(results,solution)) {
        if (verbose)
        {
            cout << "Matrix " << operation << " passed" << endl;
        }
        return true;
    }
    else {
        cout << "Matrix " << operation << " failed" << endl;
        if (verbose)
        {
            cout << "Results = " << endl;
            results.print();
            cout << "Solution = " << endl;
            solution.print();
        }
        return false;
    }
}

bool TestSuite::operationSU3Test(bool verbose, SU3 results, SU3 solution, std::string operation)
{
    if (compareSU3(results,solution)) {
        if (verbose)
        {
            cout << "Matrix " << operation << " passed" << endl;
        }
        return true;
    }
    else {
        cout << "Matrix " << operation << " failed" << endl;
        if (verbose)
        {
            cout << "Results = " << endl;
            results.print();
            cout << "Solution = " << endl;
            solution.print();
        }
        return false;
    }
}

// 2x2 MATRIX TESTS =======================================
bool TestSuite::testSU2Addition(bool verbose)
{
    return operationSU2Test(verbose,s1+s2,sAdd,"addition");
}

bool TestSuite::testSU2Subtraction(bool verbose)
{
    return operationSU2Test(verbose,s1-s2,sSub,"subtraction");
}

bool TestSuite::testSU2Multiplication(bool verbose)
{
    return operationSU2Test(verbose,s1*s2,sMul,"multiplication");
}

bool TestSuite::testSU2Transpose(bool verbose)
{
    s3.copy(s1);
    return operationSU2Test(verbose,s3.transpose(),sT,"transpose");
}

bool TestSuite::testSU2Conjugation(bool verbose)
{
    s3.copy(s1);
    return operationSU2Test(verbose,s3.conjugate(),sC,"conjugate");
}

bool TestSuite::testSU2ComplexConjugation(bool verbose)
{
    s3.copy(s1);
    s3.conjugate();
    return operationSU2Test(verbose,s3.transpose(),sCT,"conjugate transpose");
}

// 3x3 MATRIX TESTS =======================================
bool TestSuite::testSU3Addition(bool verbose)
{
    return operationSU3Test(verbose,U1+U2,UAdd,"addition");
}

bool TestSuite::testSU3Subtraction(bool verbose)
{
    return operationSU3Test(verbose,U1-U2,USub,"subtraction");
}

bool TestSuite::testSU3Multiplication(bool verbose)
{
    return operationSU3Test(verbose,U1*U2,UMul,"multiplication");
}

bool TestSuite::testSU3Transpose(bool verbose)
{
    U3.copy(U1);
    return operationSU3Test(verbose,U3.transpose(),UT,"transpose");
}

bool TestSuite::testSU3Conjugation(bool verbose)
{
    U3.copy(U1);
    return operationSU3Test(verbose,U3.conjugate(),UC,"conjugate");
}

bool TestSuite::testSU3ComplexConjugation(bool verbose)
{
    U3.copy(U1);
    U3.conjugate();
    return operationSU3Test(verbose,U3.transpose(),UCT,"conjugate transpose");
}

// SU2 PROPERTIES TESTS ===================================
// Create sub-functions that can take a SU3/SU2 matrix as an parameter and test that manually.
bool TestSuite::testSU2Hermicity(bool verbose)
{

}

bool TestSuite::testSU2Orthogonality(bool verbose)
{

}

bool TestSuite::testSU2Norm(bool verbose)
{

}

bool TestSuite::testSU2Determinant(bool verbose)
{

}

// SU3 PROPERTIES TESTS ===================================
bool TestSuite::testSU3Hermicity(bool verbose)
{
    bool passed = true;
    if (!checkSU3Hermicity(verbose, m_SU3Generator->generateRandom())) {
        cout << "FAILED: SU3 random generator matrix is not hermitian." << endl;
        passed = false;
    }
    if (!checkSU3Hermicity(verbose, m_SU3Generator->generateRST())) {
        cout << "FAILED: SU3 RST generator matrix is not hermitian." << endl;
        passed = false;
    }
    if (passed) {
        cout << "SUCCESS: SU3 random matrices are hermitian." << endl;
    }
    return passed;
}

bool TestSuite::checkSU3Hermicity(bool verbose, SU3 H)
{
    /*
     * Checks the hermicity of the SU3 matrices.
     */
    bool passed = true;
    double eps = 1e-14;
    SU3 I;
    I = H*H.inv();
    for (int i = 0; i < 18; i++) {
        if (i==0 || i==8 || i==16) continue;
        if (fabs(I[i]) > eps) passed = false;
    }
    if ((fabs(I[0] - 1) > eps) || (fabs(I[8] - 1) > eps) || (fabs(I[16] - 1) > eps)) passed = false;

    if (passed) {
        return true;
    }
    else {
        if (verbose) {
            cout << "H = " << endl;
            H.print();
            cout << "H^-1 = " << endl;
            H.inv().print();
            cout << "H^-1*H = " << endl;
            I.print();
        }
        return false;
    }
}

bool TestSuite::testSU3Orthogonality(bool verbose)
{

}

bool TestSuite::testSU3Norm(bool verbose)
{
//    double sumVal = 0; // TESTS IF WE HAVE THE NORM = 1
//    for (int i = 0; i < 3; i++) {
//        sumVal += H.normSquared(6*i);
//    }
//    cout << sqrt(sumVal) << endl;
//    exit(1);
}

bool TestSuite::testSU3Determinant(bool verbose)
{

}

// Other tests
bool TestSuite::testSU3TraceMultiplication(bool verbose)
{

}

bool TestSuite::testRSTMultiplication(bool verbose)
{

}

