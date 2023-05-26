#include "su3operations.h"
#include "parallelization/communicator.h"
#include <iomanip>

SU3Operations::SU3Operations()
{

}

// Matrix operation tester
bool SU3Operations::operationSU3Test(SU3 results, SU3 solution, std::string operation)
{
    /*
     * Function that checks if a opertion on a 3x3 matrix has succeded.
     */
    if (compareSU3(results,solution)) {
        if (m_verbose)
        {
            cout << "    SUCCESS: Matrix " << operation << " passed" << endl;
        }
        return true;
    }
    else {
        cout << "    FAILED: Matrix " << operation << " failed" << endl;
        if (m_verbose)
        {
            cout << "Results = " << endl;
            results.print();
            cout << "Solution = " << endl;
            solution.print();
        }
        return false;
    }
}
////////////////////////////////////
///////// 3x3 MATRIX TESTS /////////
////////////////////////////////////
bool SU3Operations::testSU3Addition()
{
    return operationSU3Test(U1+U2,UAdd,"addition");
}

bool SU3Operations::testSU3Subtraction()
{
    return operationSU3Test(U1-U2,USub,"subtraction");
}

bool SU3Operations::testSU3Multiplication()
{
    return operationSU3Test(U1*U2,UMul,"multiplication");
}

bool SU3Operations::testSU3Transpose()
{
    U3 = U1;
    return operationSU3Test(U3.transpose(),UT,"transpose");
}

bool SU3Operations::testSU3Conjugation()
{
    U3 = U1;
    return operationSU3Test(U3.conjugate(),UC,"conjugate");
}

bool SU3Operations::testSU3ComplexConjugation()
{
    U3 = U1;
    U3.conjugate();
    return operationSU3Test(U3.transpose(),UCT,"conjugate transpose");
}

////////////////////////////////////
/////// SU3 PROPERTIES TESTS ///////
////////////////////////////////////
bool SU3Operations::testSU3Hermicity()
{
    /*
     * Overarching function for checking the hermicity of the SU2 matrix.
     */
    bool passed = true;
    if (!checkSU3Hermicity(m_SU3Generator->generateRandom())) {
        cout << "    FAILED: SU3 random generator matrix is not hermitian." << endl;
        passed = false;
    }
    if (!checkSU3Hermicity(m_SU3Generator->generateRST())) {
        cout << "    FAILED: SU3 RST generator matrix is not hermitian." << endl;
        passed = false;
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: random SU3 matrices generated with RST and random methods are hermitian." << endl;
    }
    return passed;
}

bool SU3Operations::checkSU3Hermicity(SU3 H)
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
        if (m_verbose) {
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

bool SU3Operations::testSU3Orthogonality()
{
    /*
     * Checks if SU3 matrices generated randomly by RST or Random is orthogonal.
     */
    bool passed = true;
    if (!checkSU3Orthogonality(m_SU3Generator->generateRandom())) {
        cout << "    FAILED: SU3 random generator matrix is not orthogonal." << endl;
        passed = false;
    }
    if (!checkSU3Orthogonality(m_SU3Generator->generateRST())) {
        cout << "    FAILED: SU3 RST generator matrix is not orthogonal." << endl;
        passed = false;
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: random SU3 matrices generated with RST and random methods are orthogonal matrices." << endl;
    }
    return passed;
}

bool SU3Operations::checkSU3Orthogonality(SU3 H)
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
    if (m_verbose && !passed) {
        cout << "Column 1 and 2: " << std::setprecision(10) << c12dot << endl;
        cout << "Column 1 and 3: " << std::setprecision(10) << c13dot << endl;
        cout << "Column 2 and 3: " << std::setprecision(10) << c23dot << endl;
    }
    return passed;
}

bool SU3Operations::testSU3Norm()
{
    /*
     * Overarching function for testing if the columns of the SU3 matrices
     * generated by the RST and Random method are at unity.
     */
    bool passed = true;
    if (!checkSU3Norm(m_SU3Generator->generateRandom())) {
        cout << "    FAILED: columns of SU3 random generator matrix are not of length 1." << endl;
        passed = false;
    }
    if (!checkSU3Norm(m_SU3Generator->generateRST())) {
        cout << "    FAILED: columns of SU3 RST generator matrix are not of length 1." << endl;
        passed = false;
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: random SU3 matrices generated with RST and random methods have columns of length 1." << endl;
    }
    return passed;
}

bool SU3Operations::checkSU3Norm(SU3 H)
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
            if (m_verbose) cout << "    FAILED: norm is not zero: " << sum << endl;
            return false;
        }
    }
    return passed;
}

bool SU3Operations::testSU3Determinant()
{
    /*
     * Overarching function for testing the SU3 determinant.
    */
    bool passed = true;
    if (!checkSU3Determinant(m_SU3Generator->generateRandom())) {
        cout << "    FAILED: determinant of SU3 random generator matrix is not 1." << endl;
        passed = false;
    }
    if (!checkSU3Determinant(m_SU3Generator->generateRST())) {
        cout << "    FAILED: determinant of SU3 RST generator matrix is not 1." << endl;
        passed = false;
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: random SU3 matrix determinants is 1." << endl;
    }
    return passed;
}

bool SU3Operations::checkSU3Determinant(SU3 H)
{
    /*
     * Function that checks the determinant of the SU2 matrix.
     */
    bool passed = true;
    complex det = SU3Determinant(H);
    if (!((fabs(det.re()) - 1) < m_eps) && !(det.im() < m_eps)) {
        passed = false;
        if (m_verbose) cout << "    FAILED: the determinant of the SU3 matrix differs from 1: " << det << endl;
    }
    return passed;
}

bool SU3Operations::testSU3Trace()
{
    /*
     * Tests the trace function of the SU3 matrix class.
     */
    bool passed = true;
    complex U1TraceResult(6,6);
    if (!compareComplex(U1.trace(),U1TraceResult)) {
        passed = false;
        cout << "    FAILED: trace() of a SU3 matrix is not correct." << endl;
    }
    if (m_verbose && passed) cout << "    SUCCESS: SU3 trace is correct." << endl;
    return passed;
}

bool SU3Operations::testHermitian()
{
    /*
     * Tests that a matrix becomes hermitian if .makeHermitian() is used.
     * Translates an anti-hermitian matrix to a hermitian one.
     */
    bool passed = true;
    SU3 U = U1AntiHermitian;
    U.makeHermitian();
    for (int i = 0; i < 18; i++) {
        if (fabs(U.mat[i] - U1Hermitian.mat[i]) > 1e-15) {
            passed = false;
            cout << "    FAILED: makeHermitian() of a SU3 matrix is not correct." << endl;
            break;
        }
    }
    if (m_verbose && passed) cout << "    SUCCESS: makeHermitian() is correct." << endl;
    return passed;
}

bool SU3Operations::testAntiHermitian()
{
    /*
     * Tests that a matrix becomes anti-hermitian if .makeAntiHermitian() is used.
     * Translates an hermitian matrix to a anti-hermitian one.
     */
    bool passed = true;
    SU3 U = U1Hermitian;
    U.makeAntiHermitian();
    for (int i = 0; i < 18; i++) {
        if (fabs(U.mat[i] - U1AntiHermitian.mat[i]) > 1e-15) {
            passed = false;
            cout << "    FAILED: makeAntiHermitian() of a SU3 matrix is not correct." << endl;
            break;
        }
    }
    if (m_verbose && passed) cout << "    SUCCESS: makeAntiHermitian() is correct." << endl;
    return passed;
}

////////////////////////////////////
/////////// OTHER TESTS ////////////
////////////////////////////////////
bool SU3Operations::testSU3TraceMultiplication()
{
    /*
     * Function for ensuring that the trace multiplication performed in the WilsonGaugeAction class is correct.
     * Results retrieved from performing a simple calculation with numpy in python.
     */
    bool passed = true;
    double results = traceRealMultiplication(U1,UTrace);
    if (fabs(results - m_tracedMatrix) < m_machineEpsilon) {
        if (m_verbose) cout << "    SUCCESS: traced real results of matrix multiplication are correct." << endl;

    } else {
        if (m_verbose) cout << results << endl;
        cout << "    FAILED: traced real results of matrix multiplication are wrong." << endl;
        passed = false;
    }
    return passed;
}

// Test 3 other trace functions as well

bool SU3Operations::run3x3MatrixTests()
{
    /*
     * Function that runs tests on folowing properties of complex 2x2 matrices:
     *  - addition
     *  - subtraction
     *  - multiplication
     *  - conjugation
     *  - matrix transpose
     */
    if (m_verbose) cout << "Running basic SU3 matrix tests." << endl;
    bool passed = (testSU3Addition() && testSU3Subtraction() && testSU3Multiplication()
                   && testSU3Conjugation() && testSU3Transpose());
    if (passed) {
        cout << "PASSED: Basic SU3 matrix operations." << endl;
    } else {
        cout << "FAILURE: Basic SU3 matrix operations." << endl;
    }
    return passed;
}


bool SU3Operations::runSU3Tests()
{
    /*
     * Function that runs tests on folowing properties of SU2 matrices:
     *  - hermicity
     *  - orthogonality of columns
     *  - normality of columns
     *  - determinant equals abs(1)
     */
    bool passed = false;

    if (m_processRank == 0) {
        passed = run3x3MatrixTests();

        if (m_verbose) cout << "Running SU3 property tests." << endl;

        passed = (passed && testSU3Hermicity() && testSU3Orthogonality()
                  && testSU3Norm() && testSU3Determinant() && testAntiHermitian()
                  && testHermitian() && testSU3Trace() && testSU3TraceMultiplication());

        if (passed) {
            cout << "PASSED: SU3 properties." << endl;
        } else {
            cout << "FAILURE: SU3 properties." << endl;
        }
    }

    MPI_Bcast(&passed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    Parallel::Communicator::setBarrier();

    return passed;
}
