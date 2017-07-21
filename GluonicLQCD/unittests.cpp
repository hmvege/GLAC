#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "matrices/su3.h"
#include "matrices/su3matrixgenerator.h"
#include "correlators/plaquette.h"
#include "links.h"
#include "complex.h"
#include "functions.h"

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

void SU3BaseTests()
{
    /*
     * Basic tests of the complex class.
     */
    complex a = complex(1,1);
    complex b = complex(2,2);
    a *= b;
    cout << a << endl;
    cout << a.conjugate() << endl;
    SU3 A, B, C, D;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A.mat[(i*3+j)] = complex(i,j);
            B.mat[(i*3+j)] = complex(i*i,j*j);
            C.mat[(i*3+j)] = complex(0,0);
            D.mat[(i*3+j)] = complex(0,0);
        }
    }

    cout << endl;
    A.print();
    cout << endl;
    B.print();
    cout << A.get(1,1) << endl;
    C = A + B;
    D = A * B;
    D.print();
    cout << endl;
    cout << "Printing the first column of matrix D: " << endl;
    for (int i = 0; i < 3; i++) {
        cout << "D[" <<i<<"] = "<< D.mat[3*i].re << " + i" << D.mat[i].im << endl;
    }
    cout << endl;
    cout << "SU3 base test completed." << endl;
}

complex dot(complex * a, complex * b) {
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

bool testOrthogonality(SU3 H, bool verbose)
{
    /*
     * Small test for testing the orthogonatility of a SU3 matrix M.
     * H =
     * 0 1 2
     * 3 4 5
     * 6 7 8
     */
    complex *col1 = new complex[3];
    complex *col2 = new complex[3];
    complex *col3 = new complex[3];
    for (int i = 0; i < 3; i++) {
        col1[i] = H[3*i];
        col2[i] = H[3*i+1];
        col3[i] = H[3*i+2];
    }

    double eps = 1e-15;
    bool testPassed = true;
    complex c12dot = dot(col1,col2);
    complex c13dot = dot(col1,col3);
    complex c23dot = dot(col2,col3);
    if (verbose) {
        cout << "Column 1 and 2: " << c12dot << endl;
        cout << "Column 1 and 3: " << c13dot << endl;
        cout << "Column 2 and 3: " << c23dot << endl;
    }
    if ((fabs(c12dot.re) > eps) || (fabs(c12dot.im) > eps)) testPassed = false;
    if ((fabs(c13dot.re) > eps) || (fabs(c13dot.im) > eps)) testPassed = false;
    if ((fabs(c23dot.re) > eps) || (fabs(c23dot.im) > eps)) testPassed = false;

    if (testPassed) {
        cout << "PASSED: Columns is orthogonal." << endl;
        return true;
    }
    else {
        cout << "FAILED: Columns is not orthogonal." << endl;
        return false;
    }
}

bool testHermicity(SU3 H, bool verbose)
{
    /*
     * Small test for testing if we have Hermicity.
     */
    bool testPassed = true;
    double eps = 1e-15;
    if (verbose) {
        cout << "Matrix = " << endl;
        H.print();
    }
    SU3 I;
    I = H*inverse(H);
    if (verbose) {
        cout << "\nInverse matrix = "<< endl;
        inverse(H).print();
        cout << "\nUnit matrix = " << endl;
        I.print();
        cout << endl;
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; i < j; j++) {
            if ((fabs(I[i*3+j].re) > eps) || (fabs(I[i*3+j].im) > eps)) {
                testPassed = false;
            }
            if ((fabs(I[j*3+i].re) > eps) || (fabs(I[j*3+i].im) > eps)) {
                testPassed = false;
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        if ((fabs(I[3*i+i].re - 1) > eps) || (fabs(I[3*i+i].im) > eps)) {
            testPassed = false;
        }
    }
    if (testPassed) {
//        cout << "PASSED: matrix is hermitian." << endl;
        return true;
    }
    else {
        cout << "FAILED: matrix is not hermitian." << endl;
        I.print();
        exit(1);
        return false;
    }
}

bool testNorm(int col, SU3 H)
{
    // TEST: Finding length of column
    double s=0;
    for (int i = 0; i < 3; i++) {
        s += H[3*i+col].norm();
    }
    if (fabs(s-1) < 1e-15) {
        cout << "PASSED: length = " << s << endl;
        return true;
    }
    else {
        cout << "FAILED: length = " << setprecision(16) << s << endl;
        return false;
    }
}

void testMatrixMultiplication()
{
    SU3 A,B,C;
    // Setting A values
    A.mat[0].re = 1;
    A.mat[1].re = 0;
    A.mat[2].re = 0;
    A.mat[3].re = 0;
    A.mat[4].re = 0;
    A.mat[5].re = 0;
    A.mat[6].re = 1;
    A.mat[7].re = 0;
    A.mat[8].re = 0;
    // Setting B values
    B.mat[0].re = 0;
    B.mat[1].re = 1;
    B.mat[2].re = 0;
    B.mat[3].re = 0;
    B.mat[4].re = 0;
    B.mat[5].re = 0;
    B.mat[6].re = 0;
    B.mat[7].re = 0;
    B.mat[8].re = 0;
    C = A*B;
    printf("Matrix A: \n");
    A.print();
    printf("Matrix B: \n");
    B.print();
    printf("Matrix C: \n");
    C.print();
}

void testMatrixSU3Properties() {
    /* Testing the basic SU3 matrix properties for a known result:
     * Original:
     * U1 =
     * [[ 1.+1.j  1.+2.j  1.+3.j]
     *  [ 2.+1.j  2.+2.j  2.+3.j]
     *  [ 3.+1.j  3.+2.j  3.+3.j]]
     * U2 =
     * [[ 1.+1.j  1.+2.j  1.+3.j]
     *  [ 2.+1.j  2.+2.j  2.+3.j]
     *  [ 3.+1.j  3.+2.j  3.+3.j]]
     */
    SU3 U1, U2, U3;
    U1.mat[0] = complex(1,1);
    U1.mat[1] = complex(1,2);
    U1.mat[2] = complex(1,3);
    U1.mat[3] = complex(2,1);
    U1.mat[4] = complex(2,2);
    U1.mat[5] = complex(2,3);
    U1.mat[6] = complex(3,1);
    U1.mat[7] = complex(3,2);
    U1.mat[8] = complex(3,3);
    U2.mat[0] = complex(4,4);
    U2.mat[1] = complex(4,5);
    U2.mat[2] = complex(4,6);
    U2.mat[3] = complex(5,4);
    U2.mat[4] = complex(5,5);
    U2.mat[5] = complex(5,6);
    U2.mat[6] = complex(6,4);
    U2.mat[7] = complex(6,5);
    U2.mat[8] = complex(6,6);
    // Original matrices
    cout << "Original matrices\nU1 =" << endl;
    U1.print();
    cout << "U2 =" << endl;
    U2.print();
    // Addition
    cout << "Testing SU3 matrix addition" << endl;
    U3 = U1 + U2;
    U3.print();
    U3.zeros();
    // Subtraction
    cout << "Testing SU3 matrix subtraction" << endl;
    U3 = U1 - U2;
    U3.print();
    U3.zeros();
    // Multiplication
    cout << "Testing SU3 matrix multiplication" << endl;
    U3 = U1 * U2;
    U3.print();
    U3.zeros();
    // Conjugation
    cout << "Testing SU3 matrix conjugation" << endl;
    U3.copy(U1);
    U3.conjugate();
    U3.print();
    U3.zeros();
    // Transpose
    cout << "Testing SU3 matrix transposing" << endl;
    U3.copy(U1);
    U3.transpose();
    U3.zeros();
    // Conjugate transpose
    cout << "Testing SU3 matrix conjugate transpose" << endl;
    U3.copy(U1);
    U3.conjugate();
    U3.transpose();
    U3.print();
}

bool SU2UnitTest(complex * r, complex * s, complex * t)
{
    /*
     * Ensuring that the inverse is the transposed complex conjugated in order to produce the inverse.
     */
    complex * r_inv = new complex[4];
    complex * s_inv = new complex[4];
    complex * t_inv = new complex[4];
    complex * I = new complex[4];
    r_inv[0] = r[0].c();
    r_inv[1] = r[2].c();
    r_inv[2] = r[1].c();
    r_inv[3] = r[3].c();

    s_inv[0] = s[0].c();
    s_inv[1] = s[2].c();
    s_inv[2] = s[1].c();
    s_inv[3] = s[3].c();
    t_inv[0] = t[0].c();
    t_inv[1] = t[2].c();
    t_inv[2] = t[1].c();
    t_inv[3] = t[3].c();
    // Creating inverse matrix
    I[0] = r[0]*r_inv[0] + r[1]*r_inv[2];
    I[1] = r[0]*r_inv[1] + r[1]*r_inv[3];
    I[2] = r[2]*r_inv[0] + r[3]*r_inv[2];
    I[3] = r[2]*r_inv[1] + r[3]*r_inv[3];

    bool testPassed = true;
    double eps = 1e-15;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; i < j; j++) {
            if ((fabs(I[i*2+j].re) > eps) || (fabs(I[i*2+j].im) > eps)) {
                testPassed = false;
            }
            if ((fabs(I[j*2+i].re) > eps) || (fabs(I[j*2+i].im) > eps)) {
                testPassed = false;
            }
        }
    }

    for (int i = 0; i < 2; i++) {
        if ((fabs(I[2*i+i].re - 1) > eps) || (fabs(I[2*i+i].im) > eps)) {
            testPassed = false;
        }
    }
    delete [] r_inv;
    delete [] s_inv;
    delete [] t_inv;
    if (testPassed) {
        cout << "PASSED: SU2 matrix is hermitian." << endl;
        return true;
    } else {
        cout << "FAILED: SU2 matrix is not hermitian." << endl;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                cout << std::setw(12) << I[i*2+j] << " ";
            }
            cout << endl;
        }
        delete [] I;
        exit(1);
        return false;
    }
}

void testDeterminant(SU3 U) {
    double eps = 1e-16;
    complex det = SU3Determinant(U);
    if (((det.re - 1) < eps) && (det.im < eps)) {
        cout << "PASSED: the determinant of the SU3 matrix is 1." << endl;
    } else {
        cout << "FAILED: the determinant of the SU3 matrix differs from 1." << endl;
    }
}

void runMatrixPerformanceTest(double eps, double seed, int NTests, bool testMatrix, bool testComplex) {
    /*
     * For running performance tests of the matrix multiplication contained in SU3.
     */
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> uni_dist(0,1);
    if (testMatrix) {
        SU3 U1, U2;
        SU3MatrixGenerator SU3Gen(eps, seed);
        clock_t programStart, programEnd;
        programStart = clock();
        for (int i = 0; i < NTests; i++) {
            U1 = SU3Gen.generateRST();
            U2 *= U1;
        }
        programEnd = clock();
        cout << "Matrix multiplication performance test completed. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    }
    if (testComplex) {
        complex c1(1,1),c2;
        clock_t programStart, programEnd;
        programStart = clock();
        for (int i = 0; i < NTests; i++) {
            c2.re = uni_dist(gen);
            c2.im = uni_dist(gen);
            c1 *= c2;
        }
        programEnd = clock();
        cout << "Complex multiplication performance test completed. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    }

}

//bool checkGauge(Links *lattice, int N, SU3MatrixGenerator *SU3Generator) {
//    /*
//     * Checking the gauge invariance of the lattice.
//     */
//    bool testPassed = true;
//    SU3 preM;
//    SU3 M;
//    SU3 randomU;
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < 4; j++) {
//            randomU.copy(SU3Generator->generate());
//            M.copy(lattice[i].U[j]);
//            preM.copy(M);
//            M = inverse(randomU)*preM*randomU;
//            for (int k = 0; k < 9; k++) {
//                if (M.mat[k] != preM.mat[k]) {
//                    testPassed = false;
//                    cout << "Error in gauge invariance!" << endl;
//                }
//            }
//        }
//    }
//}

void checkDim(int N, int N_T) {
    if (N == N_T) {
        cout << "Error" << endl;
        exit(1);
    }
}
