#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <random>
#include "actions/action.h"
#include "correlators/correlator.h"
#include "links.h"
#include "su3matrixgenerator.h"

class Action;

class Metropolis
{
private:
    int N;
    int NCf;
    int NCor;
    int NTherm;
    int m_nUpdates = 10; // N updates before calculating the action, as that is costly
    double epsilon;
    double a; // Lattice spacing
    double L; // Lattice size
    int acceptanceCounter = 0;

    // Lattice constants
    int latticeSize;
    Links * lattice;
    SU3 updatedMatrix;

    // Variables used to perform statistics
    double * Gamma;
    double * GammaSquared;
//    double * GammaVariance;
//    double * GammaStd;
    double averagedGamma;
    double varianceGamma;
    double stdGamma;
//    double deltaE_std;
//    double deltaE;

    // Storing the action as a pointer
    Action *S = nullptr;
    double deltaS;

    // For sampling the system(lattice)
    void sampleSystem();

    // Correlator
    Correlator * m_correlator = nullptr;

    // Function for updating our system using the Metropolis algorithm
    void update();
    void updateLink(int i, int mu);

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;
public:
    Metropolis(int new_N, int new_NCf, int new_NCor, int NTherm, double new_a, double new_L, Correlator *new_correlator, Action *new_S);
    ~Metropolis();
    void latticeSetup(SU3MatrixGenerator *SU3Generator);
    void runMetropolis();
    void getStatistics();
//    void writeStatisticsToFile(const char *filename);
    void writeDataToFile(const char *filename);

    // Setters
    void setAction(Action *newS) { S = newS; }
    void setCorrelator(Correlator *correlator) { m_correlator = correlator; }
    void setN(int new_N) { N = new_N; }
    void setNCf(int new_NCf) { NCf = new_NCf; }
    void setEpsilon(double new_epsilon) { epsilon = new_epsilon; }
    void setUpdateFrequency(int nUpdates) { m_nUpdates = nUpdates; }

    // Getters
    int getN() { return N; }
    int getNCf() { return NCf; }
    int getEpsilon() { return epsilon; }

    // Printers
    void printEnergies();
    void printAcceptanceRate();
};

#endif // METROPOLIS_H
