#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "math/matrices/su3.h"
#include "math/matrices/su3matrixgenerator.h"
#include "observables/plaquette.h"
#include "math/links.h"
#include "math/complex.h"
#include "math/functions.h"

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
    complex c = complex(1,1);
    complex d = complex(2,2);
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    a *= b;
    cout << "a*b = " << a << endl;
//    cout << "conjugate(a) = " << a.c() << endl;
    c += d;
    cout << "1+1j + 2+2j = " << c << endl;

    a /= b;
    cout << "a/b = " << a << endl;

//    SU3 A, B, C, D;
//    for (int i = 0; i < 3; i++)
//    {
//        for (int j = 0; j < 3; j++)
//        {
//            A.setComplex(complex(i,j), 3*i+j);
//            B.setComplex(complex(i*i,j*j), 3*i+j);
//            C.setComplex(complex(0,0), 3*i+j);
//            D.setComplex(complex(0,0), 3*i+j);
////            A.mat[(i*3+j)] = complex(i,j);
////            B.mat[(i*3+j)] = complex(i*i,j*j);
////            C.mat[(i*3+j)] = complex(0,0);
////            D.mat[(i*3+j)] = complex(0,0);
//        }
//    }

//    cout << endl;
//    cout << "Printing A" << endl;
//    A.print();
//    cout << endl;
//    cout << "Printing B" << endl;
//    B.print();
//    cout << "Getting middle element of A" << endl;
//    cout << A.get(1,1) << endl;
//    C += A;
//    C = A + B;
//    D = A * B;
//    cout << "Printing D" << endl;
//    D.print();
//    cout << endl;
//    cout << "Printing the first column of matrix D: " << endl;
//    for (int i = 0; i < 3; i++) {
//        cout << "D[" <<i<<"] = "<< D.mat[6*i] << " + i" << D.mat[6*i+1] << endl;
//    }
//    cout << endl;
//    cout << "SU3 base test completed." << endl;
}


void runMatrixPerformanceTest(double seed, int NTests, bool testMatrix, bool testComplex) {
    /*
     * For running performance tests of the matrix multiplication contained in SU3.
     */
    cout << "Starting matrix performance test." << endl;
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> uni_dist(0,1);
    if (testMatrix) {
        SU3 U1, U2;
        SU3MatrixGenerator SU3Gen;
        clock_t programStart, programEnd;
        programStart = clock();
        for (int i = 0; i < NTests; i++) {
            U1 = SU3Gen.generateRST();
            U2 = SU3Gen.generateRST();
            U2 *= U1;
//            U1 = SU3Gen.generateRandom();
//            U1.print();
//            U2 = U1.antiHermitian();
//            U2.print();
//            cout << U2.trace() << endl;
//            (U1*U2).print();
//            U2 *= U1;
        }
        programEnd = clock();
        cout << "Matrix multiplication performance test completed. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    }
    if (testComplex) {
        complex c1(1,1),c2;
        clock_t programStart, programEnd;
        programStart = clock();
        for (int i = 0; i < NTests; i++) {
            c2.setRe(uni_dist(gen));
            c2.setIm(uni_dist(gen));
            c1 *= c2;
        }
        programEnd = clock();
        cout << "Complex multiplication performance test completed. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    }
}

void runBoolTest(int NTests) {
    /*
     * Small tester for certain functions performance.
     */

    std::mt19937_64 gen(-1);
    std::uniform_int_distribution<int> int_dist(0,1);
    int int1 = 0, int2 = 0;
    clock_t programStart, programEnd;
    long unsigned int counter = 0;
    programStart = clock();
    for (int i = 0; i < NTests; i++)  {
        int1 = int_dist(gen);
        int2 = int_dist(gen);
        if (int1 && int2) {
            counter++;
        }
    }
    programEnd = clock();
    cout << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
}
