#ifndef SU3MATRIXGENERATOR_H
#define SU3MATRIXGENERATOR_H

#include <random>
#include "su3.h"
#include "su2.h"
#include "complex.h"

class SU3MatrixGenerator
{
private:
    double epsilon;
    double epsilonSquared;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> uniform_distribution;
    std::uniform_real_distribution<double> SU2_uniform_distribution;
    void GramSchmitt();

    // Used for creating a random matrix
    SU3 H;
    // Used for creating a random matrix close to unity
    SU3 X, R, S, T;
    SU2 r,s,t;

    // Used for generating SU2 matrices
    SU2 U;
    double _r[4];
    double _x[4];
    double _rNorm = 0;

    SU2 *sigma;
public:
    SU3MatrixGenerator(double eps, double seed);
    ~SU3MatrixGenerator();
    SU3 generateRandom();
    SU3 generateRST();
    SU2 generateSU2();

    // Setters
    void setEpsilon(double eps) { epsilon = eps; }

    // Getters
    double getEpsilon() { return epsilon; }
};

#endif // SU3MATRIXGENERATOR_H
