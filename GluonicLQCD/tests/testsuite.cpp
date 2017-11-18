#include "testsuite.h"

using std::cout;
using std::endl;

TestSuite::TestSuite()
{
    /*
     * Class for running unit tests.
     */
    // Initiating the Mersenne Twister random number generator
    m_generator = std::mt19937_64(std::time(nullptr));
    m_uniform_distribution = std::uniform_real_distribution<double>(0,1); // Print out max values to ensure we dont go out of scope!!
    // Initiating the SU3 Matrix generator
    m_SU3Generator = new SU3MatrixGenerator;

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
    // RST Matrix multiplier
    U_RST.setComplex(complex(-8,9),0);
    U_RST.setComplex(complex(-123,82),2);
    U_RST.setComplex(complex(-78,83),4);
    U_RST.setComplex(complex(-20,37),6);
    U_RST.setComplex(complex(-391,342),8);
    U_RST.setComplex(complex(-234,319),10);
    U_RST.setComplex(complex(2,6),12);
    U_RST.setComplex(complex(48,58),14);
    U_RST.setComplex(complex(44,35),16);
    // r SU2 matrix in RST multiplier test
    s_r.setComplex(complex(1,2),0);
    s_r.setComplex(complex(3,4),2);
    s_r.setComplex(complex(5,6),4);
    s_r.setComplex(complex(7,8),6);
    // s SU2 matrix in RST multiplier test
    s_s.setComplex(complex(2,5),0);
    s_s.setComplex(complex(4,8),2);
    s_s.setComplex(complex(2,6),4);
    s_s.setComplex(complex(10,3),6);
    // t SU2 matrix in RST multiplier test
    s_t.setComplex(complex(7,2),0);
    s_t.setComplex(complex(6,1),2);
    s_t.setComplex(complex(6,4),4);
    s_t.setComplex(complex(5,2),6);
    // Traced matrix multiplication
    m_tracedMatrix = -21.0; // re(tr(U1*U3))
    UTrace.setComplex(complex(4,4),0);
    UTrace.setComplex(complex(5,5),2);
    UTrace.setComplex(complex(2,6),4);
    UTrace.setComplex(complex(6,4),6);
    UTrace.setComplex(complex(5,5),8);
    UTrace.setComplex(complex(1,6),10);
    UTrace.setComplex(complex(6,4),12);
    UTrace.setComplex(complex(9,5),14);
    UTrace.setComplex(complex(2,6),16);
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
    /*
     * Function that runs all available tests.
     */
    for (int i = 0; i < 60; i++) cout << "=";
    cout << endl;
    if (verbose)
    {
        // Original matrices
        cout << "Original matrices\nU1 =" << endl;
        U1.print();
        cout << "U2 =" << endl;
        U2.print();
    }
    if (run3x3MatrixTests(verbose) & run2x2MatrixTests(verbose) & runSU2Tests(verbose) & runSU3Tests(verbose) & runFunctionsTest(verbose)) {
        cout << "SUCCESS: All tests passed." << endl;
    } else {
        cout << "FAILURE: One or more test(s) failed." << endl;
    }
    for (int i = 0; i < 60; i++) cout << "=";
    cout << endl;
}

bool TestSuite::runSU2Tests(bool verbose)
{
    /*
     * Function that runs tests on folowing properties of SU2 matrices:
     *  - hermicity
     *  - orthogonality of columns
     *  - normality of columns
     *  - determinant equals abs(1)
     */
    if (verbose) cout << "Running SU2 property tests." << endl;
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
    /*
     * Function that runs tests on folowing properties of SU2 matrices:
     *  - hermicity
     *  - orthogonality of columns
     *  - normality of columns
     *  - determinant equals abs(1)
     */
    if (verbose) cout << "Running SU3 property tests." << endl;
    bool passed = (testSU3Hermicity(verbose) && testSU3Orthogonality(verbose) && testSU3Norm(verbose)
                   && testSU3Determinant(verbose));
    if (passed) {
        cout << "PASSED: SU3 properties." << endl;
    } else {
        cout << "FAILED: SU3 properties." << endl;
    }

    return passed;
}

bool TestSuite::run2x2MatrixTests(bool verbose)
{
    /*
     * Function that runs tests on folowing properties of complex 2x2 matrices:
     *  - addition
     *  - subtraction
     *  - multiplication
     *  - conjugation
     *  - matrix transpose
     */
    if (verbose) cout << "Running basic SU2 matrix tests." << endl;
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
    /*
     * Function that runs tests on folowing properties of complex 2x2 matrices:
     *  - addition
     *  - subtraction
     *  - multiplication
     *  - conjugation
     *  - matrix transpose
     */
    if (verbose) cout << "Running basic SU3 matrix tests." << endl;
    bool passed = (testSU3Addition(verbose) && testSU3Subtraction(verbose) && testSU3Multiplication(verbose)
                   && testSU3Conjugation(verbose) && testSU3Transpose(verbose));
    if (passed) {
        cout << "PASSED: Basic SU3 matrix operations." << endl;
    } else {
        cout << "FAILED: Basic SU3 matrix operations." << endl;
    }
    return passed;
}

bool TestSuite::runFunctionsTest(bool verbose)
{
    /*
     * Function that tests key elements of the program, that is:
     *  - trace multiplication of the getDeltaAction in WilsonGaugeAction class.
     *  - RST multiplication in the generateRST in SU3MatrixGenerator class.
     */
    if (verbose) cout << "Running test on specific matrix trace multiplication." << endl;
    bool passed = true;
    if (testSU3TraceMultiplication(verbose)) {
        cout << "PASSED: matrix trace multiplication test passed." << endl;
    } else {
        cout << "FAILED: matrix trace multiplication test failed." << endl;
        passed = false;
    }
    if (testRSTMultiplication(verbose)) {
        cout << "PASSED: matrix RST multiplication test passed." << endl;
    } else {
        cout << "FAILED: matrix RST multiplication test failed." << endl;
        passed = false;
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

// Inline dot product between two columns of a 3x3 matrix
inline complex TestSuite::dot(complex * a, complex * b) {
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex returnSum(0,0);
    for (int i = 0; i < 3; i++) {
        returnSum += a[i].c()*b[i];
    }
    return returnSum;
}

// Inline dot product between two columns of a 2x2 matrix
inline complex TestSuite::dot2(complex * a, complex * b) {
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex returnSum(0,0);
    for (int i = 0; i < 2; i++) {
        returnSum += a[i].c()*b[i];
    }
    return returnSum;
}

// Matrix operation tester ================================
bool TestSuite::operationSU2Test(bool verbose, SU2 results, SU2 solution, std::string operation)
{
    /*
     * Function that checks if a opertion on a 2x2 matrix has succeded.
     */
    if (compareSU2(results,solution)) {
        if (verbose)
        {
            cout << "SUCCESS: Matrix " << operation << " passed" << endl;
        }
        return true;
    }
    else {
        cout << "FAILURE: Matrix " << operation << " failed" << endl;
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
    /*
     * Function that checks if a opertion on a 3x3 matrix has succeded.
     */
    if (compareSU3(results,solution)) {
        if (verbose)
        {
            cout << "SUCCESS: Matrix " << operation << " passed" << endl;
        }
        return true;
    }
    else {
        cout << "FAILURE: Matrix " << operation << " failed" << endl;
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
bool TestSuite::testSU2Hermicity(bool verbose)
{
    /*
     * Overarching function for checking the hermicity of the SU2 matrix.
     */
    bool passed = true;
    if (!checkSU2Hermicity(verbose, m_SU3Generator->generateSU2())) {
        cout << "FAILED: randomly generated SU2 matrix is not hermitian." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: randomly generated SU2 matrix is hermitian." << endl;
    }
    return passed;
}

bool TestSuite::checkSU2Hermicity(bool verbose, SU2 H)
{
    /*
     * Checks the hermicity of the SU2 matrices.
     * H = 0 1  2 3
     *     4 5  6 7
     */
    bool passed = true;
    SU2 I;
    I = H*H.inv();
    for (int i = 0; i < 8; i++) {
        if (i==0 || i==6) continue;
        if (fabs(I[i]) > m_eps) passed = false;
    }
    if ((fabs(I[0] - 1) > m_eps) || (fabs(I[6] - 1) > m_eps)) passed = false;

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

bool TestSuite::testSU2Orthogonality(bool verbose)
{
    /*
     * Checks if SU2 matrices generated randomly by RST or Random is orthogonal.
     */
    bool passed = true;
    if (!checkSU2Orthogonality(verbose, m_SU3Generator->generateSU2())) {
        cout << "FAILED: SU2 randomly generated matrix is not orthogonal." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: randomly generated SU2 matrix is orthogonal." << endl;
    }
    return passed;
}

bool TestSuite::checkSU2Orthogonality(bool verbose, SU2 H)
{
    /*
     * Testing the orthogonatility of a SU3 matrix H.
     */
    bool passed = true;
    complex col1[2];
    complex col2[2];
    for (int i = 0; i < 2; i++) {
        col1[i].setRe(H[4*i]);
        col1[i].setIm(H[4*i+1]);
        col2[i].setRe(H[4*i+2]);
        col2[i].setIm(H[4*i+3]);
    }
    complex c12dot = dot2(col1,col2);
    if ((fabs(c12dot.re()) > m_eps) || (fabs(c12dot.im()) > m_eps)) passed = false;
    if (verbose && !passed) {
        cout << "Column 1 and 2: " << std::setprecision(10) << c12dot << endl;
    }
    return passed;
}

bool TestSuite::testSU2Norm(bool verbose)
{
    /*
     * Overarching function for testing if the columns of the SU3 matrices
     * generated by the RST and Random method are at unity.
     */
    bool passed = true;
    if (!checkSU2Norm(verbose, m_SU3Generator->generateSU2())) {
        cout << "FAILED: columns of a randomly generated SU2 matrix is not at length 1." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: randomly generated SU2 matrices have columns of length 1." << endl;
    }
    return passed;
}

bool TestSuite::checkSU2Norm(bool verbose, SU2 H)
{
    /*
     * Function that checks the SU3 norm of a matrix H.
     */
    bool passed = true;
    double sum;
    for (int col = 0; col < 2; col++) {
        sum = 0;
        for (int i = 0; i < 2; i++) {
            sum += H.normSquared(4*i + 2*col);
        }
        sum = sqrt(sum);
        if (fabs(sqrt(sum) - 1) > m_eps) {
            if (verbose) cout << "FAILED: norm is not zero: " << sum << endl;
            return false;
        }
    }
    return passed;
}

bool TestSuite::testSU2Determinant(bool verbose)
{
    /*
     * Overarching function for testing the SU2 determinant.
    */
    bool passed = true;
    if (!checkSU2Determinant(verbose, m_SU3Generator->generateSU2())) {
        cout << "FAILED: determinant of randomly generated SU2 matrix is not 1." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: random SU2 matrix determinants is 1." << endl;
    }
    return passed;
}

bool TestSuite::checkSU2Determinant(bool verbose, SU2 H)
{
    /*
     * Function that checks the determinant of the SU2 matrix.
     */
    bool passed = true;
    complex det = SU2Determinant(H);
    if (!((fabs(det.re()) - 1) < m_eps) && !(det.im() < m_eps)) {
        passed = false;
        if (verbose) cout << "FAILED: the determinant of the SU2 matrix differs from 1: " << det << endl;
    }
    return passed;
}

// SU3 PROPERTIES TESTS ===================================
bool TestSuite::testSU3Hermicity(bool verbose)
{
    /*
     * Overarching function for checking the hermicity of the SU2 matrix.
     */
    bool passed = true;
    if (!checkSU3Hermicity(verbose, m_SU3Generator->generateRandom())) {
        cout << "FAILED: SU3 random generator matrix is not hermitian." << endl;
        passed = false;
    }
    if (!checkSU3Hermicity(verbose, m_SU3Generator->generateRST())) {
        cout << "FAILED: SU3 RST generator matrix is not hermitian." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: random SU3 matrices generated with RST and random methods are hermitian." << endl;
    }
    return passed;
}

bool TestSuite::checkSU3Hermicity(bool verbose, SU3 H)
{
    /*
     * Checks the hermicity of the SU3 matrices.
     */
    bool passed = true;
    SU3 I;
    I = H*H.inv();
    for (int i = 0; i < 18; i++) {
        if (i==0 || i==8 || i==16) continue;
        if (fabs(I[i]) > m_eps) passed = false;
    }
    if ((fabs(I[0] - 1) > m_eps) || (fabs(I[8] - 1) > m_eps) || (fabs(I[16] - 1) > m_eps)) passed = false;
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
    /*
     * Checks if SU3 matrices generated randomly by RST or Random is orthogonal.
     */
    bool passed = true;
    if (!checkSU3Orthogonality(verbose, m_SU3Generator->generateRandom())) {
        cout << "FAILED: SU3 random generator matrix is not orthogonal." << endl;
        passed = false;
    }
    if (!checkSU3Orthogonality(verbose, m_SU3Generator->generateRST())) {
        cout << "FAILED: SU3 RST generator matrix is not orthogonal." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: random SU3 matrices generated with RST and random methods are orthogonal matrices." << endl;
    }
    return passed;
}

bool TestSuite::checkSU3Orthogonality(bool verbose, SU3 H)
{
    /*
     * Testing the orthogonatility of a SU3 matrix H.
     */
    bool passed = true;
    complex col1[3];
    complex col2[3];
    complex col3[3];
    for (int i = 0; i < 3; i++) {
        col1[i].setRe(H[6*i]);
        col1[i].setIm(H[6*i+1]);
        col2[i].setRe(H[6*i+2]);
        col2[i].setIm(H[6*i+3]);
        col3[i].setRe(H[6*i+4]);
        col3[i].setIm(H[6*i+5]);
    }
    complex c12dot = dot(col1,col2);
    complex c13dot = dot(col1,col3);
    complex c23dot = dot(col2,col3);
    if ((fabs(c12dot.re()) > m_eps) || (fabs(c12dot.im()) > m_eps)) passed = false;
    if ((fabs(c13dot.re()) > m_eps) || (fabs(c13dot.im()) > m_eps)) passed = false;
    if ((fabs(c23dot.re()) > m_eps) || (fabs(c23dot.im()) > m_eps)) passed = false;
    if (verbose && !passed) {
        cout << "Column 1 and 2: " << std::setprecision(10) << c12dot << endl;
        cout << "Column 1 and 3: " << std::setprecision(10) << c13dot << endl;
        cout << "Column 2 and 3: " << std::setprecision(10) << c23dot << endl;
    }
    return passed;
}

bool TestSuite::testSU3Norm(bool verbose)
{
    /*
     * Overarching function for testing if the columns of the SU3 matrices
     * generated by the RST and Random method are at unity.
     */
    bool passed = true;
    if (!checkSU3Norm(verbose, m_SU3Generator->generateRandom())) {
        cout << "FAILED: columns of SU3 random generator matrix are not of length 1." << endl;
        passed = false;
    }
    if (!checkSU3Norm(verbose, m_SU3Generator->generateRST())) {
        cout << "FAILED: columns of SU3 RST generator matrix are not of length 1." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: random SU3 matrices generated with RST and random methods have columns of length 1." << endl;
    }
    return passed;
}

bool TestSuite::checkSU3Norm(bool verbose, SU3 H)
{
    /*
     * Function that checks the SU3 norm of a matrix H.
     */
    bool passed = true;
    double sum;
    for (int col = 0; col < 3; col++) {
        sum = 0;
        for (int i = 0; i < 3; i++) {
            sum += H.normSquared(6*i + 2*col);
        }
        sum = sqrt(sum);
        if (fabs(sqrt(sum) - 1) > m_eps) {
            if (verbose) cout << "FAILED: norm is not zero: " << sum << endl;
            return false;
        }
    }
    return passed;
}

bool TestSuite::testSU3Determinant(bool verbose)
{
    /*
     * Overarching function for testing the SU3 determinant.
    */
    bool passed = true;
    if (!checkSU3Determinant(verbose, m_SU3Generator->generateRandom())) {
        cout << "FAILED: determinant of SU3 random generator matrix is not 1." << endl;
        passed = false;
    }
    if (!checkSU3Determinant(verbose, m_SU3Generator->generateRST())) {
        cout << "FAILED: determinant of SU3 RST generator matrix is not 1." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: random SU3 matrix determinants is 1." << endl;
    }
    return passed;
}

bool TestSuite::checkSU3Determinant(bool verbose, SU3 H)
{
    /*
     * Function that checks the determinant of the SU2 matrix.
     */
    bool passed = true;
    complex det = SU3Determinant(H);
    if (!((fabs(det.re()) - 1) < m_eps) && !(det.im() < m_eps)) {
        passed = false;
        if (verbose) cout << "FAILED: the determinant of the SU3 matrix differs from 1: " << det << endl;
    }
    return passed;
}

// Other tests
bool TestSuite::testSU3TraceMultiplication(bool verbose)
{
    /*
     * Function for ensuring that the trace multiplication performed in the WilsonGaugeAction class is correct.
     * Results retrieved from performing a simple calculation with numpy in python.
     */
    bool passed = true;
    double results = traceRealMultiplication(U1,UTrace);
    if (results == m_tracedMatrix) {
        if (verbose) cout << "SUCCESS: traced real results of matrix multiplication are correct." << endl;

    } else {
        if (verbose) cout << results << endl;
        cout << "FAILED: traced real results of matrix multiplication are wrong." << endl;
        passed = false;
    }
    return passed;
}

bool TestSuite::testRSTMultiplication(bool verbose)
{
    /*
     * Function for checking that we are multiplying the SU2 matrices correctly when generating a SU3 matrix close to unity.
     */
    bool passed = true;
    SU3 results = m_SU3Generator->testRSTMultiplication(s_r,s_s,s_t);
    if (compareSU3(results,U_RST)) {
        if (verbose) cout << "SUCCESS: RST multiplication test passed." << endl;

    } else {
        if (verbose) results.print();
        cout << "FAILED: RST multiplication test did not pass." << endl;
        passed = false;
    }
    return passed;
}

// Perform tests on matrix X
void TestSuite::testMatrix(SU3 X, bool verbose)
{
    bool passed = true;
    // Checks hermicity
    if (!checkSU3Hermicity(true, X)) {
        cout << "FAILED: matrix is not hermitian." << endl;
        passed = false;
    }
    // Checks orthogonality
    if (!checkSU3Orthogonality(verbose, X)) {
        cout << "FAILED: matrix is not orthogonal." << endl;
        passed = false;
    }
    // Checks the norm
    if (!checkSU3Norm(verbose, X)) {
        cout << "FAILED: columns of matrix are not of length 1." << endl;
        passed = false;
    }
    // Checks the determinant
    if (!checkSU3Determinant(verbose, X)) {
        cout << "FAILED: determinant of matrix is not 1." << endl;
        passed = false;
    }
    if (passed && verbose) {
        cout << "SUCCESS: All tests passed." << endl;
    }
}
