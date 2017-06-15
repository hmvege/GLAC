#include <iostream>
#include "su3.h"
#include "complex.h"
#include <cmath>
using std::cout;
using std::endl;
//using std::fabs;

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
    complex returnSum(0,0);
    for (int i = 0; i < 3; i++) {
        returnSum += a[i]*b[i];
    }
    return returnSum;
}

void checkOrthogonality(SU3 H)
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

    double eps = 1e-16;
    bool testPassed = true;
    complex c12dot = dot(col1,col2);
    complex c13dot = dot(col1,col3);
    complex c23dot = dot(col2,col3);
    if ((fabs(c12dot.re) > eps) || (fabs(c12dot.im) > eps)) testPassed = false;
    if ((fabs(c13dot.re) > eps) || (fabs(c13dot.im) > eps)) testPassed = false;
    if ((fabs(c23dot.re) > eps) || (fabs(c23dot.im) > eps)) testPassed = false;

    if (testPassed) {
        cout << "PASSED: Columns is orthogonal." << endl;
    }
    else {
        cout << "FAILED: Columns is not orthogonal." << endl;
    }
}
