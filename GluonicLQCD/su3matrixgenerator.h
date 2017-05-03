#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include <armadillo>

class SU3MatrixGenerator
{
private:
    int dim = 3; // Will always be 3 - not necessary?
    double epsilon;
//    arma::cx_mat complexIdentity;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_distribution;
public:
    SU3MatrixGenerator(double eps, std::mt19937_64 &gen, std::uniform_real_distribution<double> &randDistr);
    ~SU3MatrixGenerator();
    double *generate();
    void generateHermitian(arma::cx_mat &H);
    // Setters
    void setEpsilon(double eps) { epsilon = eps; }
    // Getters
    double getEpsilon() { return epsilon; }
};

#endif // SU3MATRIXGENERATOR_H
