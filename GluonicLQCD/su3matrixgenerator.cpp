#include "su3matrixgenerator.h"
#include <random>
#include <iomanip>
#include <iostream>

#include "su3.h"
#include "complex.h"

using std::cout;
using std::endl;

/*
 * Function for generating random SU3 matrices
 */

SU3MatrixGenerator::SU3MatrixGenerator(double eps, std::mt19937_64 &gen, std::uniform_real_distribution<double> &randDistr)
{
    epsilon = eps;
    generator = gen;
    uniform_distribution = randDistr;
//    complexIdentity = arma::ones<arma::cx_mat>(dim,dim);
}

SU3MatrixGenerator::~SU3MatrixGenerator()
{

}

double * SU3MatrixGenerator::generate()
{
    double * x;
//    std::cout << M(0,0) << std::endl;
    complex a = complex(1,1);
    complex b = complex(2,2);
    a *= b;
    cout << a << endl;
    SU3 A, B, C;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A.mat[(i*3+j)] = complex(i,j);
            B.mat[(i*3+j)] = complex(i*i,j*j);
            C.mat[(i*3+j)] = complex(0,0);
        }
    }
    cout << endl;
    A.print();
    cout << endl;
    B.print();
    cout << A.get(1,1) << endl;
    C = A + B;

    std::cout << "exiting in SU3MatrixGenerator.cpp line 40" << std::endl;
    exit(1);
    return x;
}

void SU3MatrixGenerator::generateHermitian()
{
    // Initializes upper triagonal part of matrix
//    for (int i = 1; i < dim; i++)
//    {
//        for (int j = i; j < dim; j++)
//        {
//            H(i,j).real(uniform_distribution(generator));
//            H(i,j).imag(uniform_distribution(generator));
//            H(j,i).real(H(i,j).real());
//            H(j,i).imag(-H(i,j).imag());
//        }
//    }
//    for (int i = 0; i < dim; i++)
//    {
//        H(i,i).real(uniform_distribution(generator));
//        H(i,i).imag(uniform_distribution(generator));
//    }
}
