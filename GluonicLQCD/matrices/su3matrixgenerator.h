#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include "su3.h"
#include "su2.h"

class SU3MatrixGenerator
{
private:
    double epsilon;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_distribution;
    std::uniform_real_distribution<double> SU2_uniform_distribution;
    void GramSchmitt();

    SU2 *sigma;
//    SU2 sigma1, sigma2, sigma3, su2Identity;
    SU2 su2Identity;
public:
    SU3MatrixGenerator(double eps, double seed);
    ~SU3MatrixGenerator();
    SU3 generate();
    SU3 generateIdentity();
    SU3 updateMatrix();
    void generateHermitian();


    SU2 generateSU2();

    // Setters
    void setEpsilon(double eps) { epsilon = eps; }

    // Getters
    double getEpsilon() { return epsilon; }
};

#endif // SU3MATRIXGENERATOR_H
