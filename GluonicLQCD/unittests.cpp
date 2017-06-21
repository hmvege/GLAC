#include <iostream>
#include <iomanip>
#include <cmath>
#include "su3.h"
#include "complex.h"

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
    SU3 HInv;
    HInv.copy(H);
    HInv.transpose();
    HInv.conjugate();
    SU3 I;
    I = H*HInv;
    if (verbose) {
        cout << "\nInverse matrix = "<< endl;
        HInv.print();
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
        cout << "PASSED: matrix is hermitian." << endl;
        return true;
    }
    else {
        cout << "FAILED: matrix is not hermitian." << endl;
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
