#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include "su3.h"

class SU3MatrixGenerator
{
private:
    double epsilon;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_distribution;
    void GramSchmitt();
public:
    SU3MatrixGenerator(double eps, double seed);
    ~SU3MatrixGenerator();
    SU3 generate();
    SU3 generateIdentity();
    SU3 updateMatrix();
    void generateHermitian();

    // Setters
    void setEpsilon(double eps) { epsilon = eps; }

    // Getters
    double getEpsilon() { return epsilon; }
};

#endif // SU3MATRIXGENERATOR_H
