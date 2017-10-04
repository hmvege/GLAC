#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include "su3.h"
#include "su2.h"

class SU3MatrixGenerator
{
private:
    double epsilon;
    double epsilonSquared;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_distribution;
    std::uniform_real_distribution<double> SU2_uniform_distribution;
    void GramSchmitt();

    // For random matrices not close to unity
    SU3 H;

    // Variables for generating SU3 matrices
    SU3 X, R, S, T;
    SU2 r,s,t;

    // Variables for generating SU2 matrices
    SU2 U;
    double _r[4];
    double _x[4];
    double _rNorm;

    SU2 *sigma;
//    SU2 sigma1, sigma2, sigma3, su2Identity;
public:
    SU3MatrixGenerator(double eps, double seed);
    ~SU3MatrixGenerator();
    SU3 generateRandom();
    SU3 generateRST();
    void generateHermitian();


    SU2 generateSU2();

    // Setters
    void setEpsilon(double eps) { epsilon = eps; }

    // Getters
    double getEpsilon() { return epsilon; }
};

#endif // SU3MATRIXGENERATOR_H
