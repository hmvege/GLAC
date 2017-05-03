#include "su3matrixgenerator.h"
#include <random>
#include <armadillo>
#include <iomanip>

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
//    arma::cx_mat H = arma::zeros<arma::mat>(dim,dim);
//    generateHermitian(H);
//    std::cout << H << std::endl;
//    arma::cx_mat unitMatrix = arma::ones<arma::mat>(dim,dim);
//    arma::cx_mat M = epsilon*H;

//    arma::mat H = arma::zeros<arma::mat>(dim,dim);
//    arma::cx_mat M = arma::cx_mat(unitMatrix,epsilon*H);
////    M(0,0).imag(10);// = 10;
    double * x;
//    std::cout << M(0,0) << std::endl;
    std::cout << "exiting in SU3MatrixGenerator" << std::endl;
    exit(1);
    return x;
}

void SU3MatrixGenerator::generateHermitian(arma::cx_mat &H)
{
    // Initializes upper triagonal part of matrix
    for (int i = 1; i < dim; i++)
    {
        for (int j = i; j < dim; j++)
        {
            H(i,j).real(uniform_distribution(generator));
            H(i,j).imag(uniform_distribution(generator));
            H(j,i).real(H(i,j).real());
            H(j,i).imag(-H(i,j).imag());
        }
    }
    for (int i = 0; i < dim; i++)
    {
        H(i,i).real(uniform_distribution(generator));
        H(i,i).imag(uniform_distribution(generator));
    }
}
