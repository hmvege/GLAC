#ifndef LINKS_H
#define LINKS_H

#include <armadillo>

struct Links
{
    // Link variables being the matrices
//    double * M0 = new double[matrixSize]; // CHANGE TO ARMADILLO
//    double * M1 = new double[matrixSize];
//    double * M2 = new double[matrixSize];
//    double * M3 = new double[matrixSize];
    double * M0; // CHANGE TO C STRUCT? Not nessecarily important, as C++ should provide similar speed
    double * M1;
    double * M2;
    double * M3;
};

#endif // LINKS_H
