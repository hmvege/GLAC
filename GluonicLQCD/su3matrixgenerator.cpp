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
}

SU3MatrixGenerator::~SU3MatrixGenerator()
{

}

void SU3MatrixGenerator::generate()
{
    double * x;

    exit(1);
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
