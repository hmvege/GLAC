#ifndef LINKS_H
#define LINKS_H

#include "su3.h"

class Links
{
public:
    Links();
    ~Links();
    // Link variables being the matrices
//    double * M0 = new double[matrixSize]; // CHANGE TO ARMADILLO
//    double * M1 = new double[matrixSize];
//    double * M2 = new double[matrixSize];
//    double * M3 = new double[matrixSize];
//    SU3 M0; // CHANGE TO C STRUCT? Not nessecarily important, as C++ should provide similar speed
//    SU3 M1;
//    SU3 M2;
//    SU3 M3;
    SU3 * U;
};

#endif // LINKS_H
