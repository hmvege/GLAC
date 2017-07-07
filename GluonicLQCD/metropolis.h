#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <random>
#include "actions/action.h"
#include "correlators/correlator.h"
#include "links.h"
#include "matrices/su3matrixgenerator.h"

class Action;

class Metropolis
{
private:
    int m_N;
    int m_N_T;
    int m_NCf;
    int m_NCor;
    int m_NTherm;
    int m_nUpdates = 10; // N updates before calculating the action, as that is costly
    double m_epsilon;
    double m_a; // Lattice spacing
    double m_L; // Lattice size
    // For handling the acceptance rate
    int acceptanceCounter = 0;
    double getAcceptanceRate();

    // Lattice constants
    int m_latticeSize;
    Links * m_lattice;
    SU3 m_updatedMatrix;

    // Variables used to perform statistics
    double * Gamma;
    double * GammaSquared;
//    double * GammaVariance;
//    double * GammaStd;
    double averagedGamma = 0;
    double varianceGamma = 0;
    double stdGamma = 0;
//    double deltaE_std;
//    double deltaE;

    // Storing the action as a pointer
    Action *m_S = nullptr;
    double m_deltaS;

    // For sampling the system(lattice)
    void sampleSystem();

    // Correlator
    Correlator * m_correlator = nullptr;

    // Function for updating our system using the Metropolis algorithm
    void update();
    void updateLink(int latticeIndex, int mu);

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;
public:
    Metropolis(int N, int N_T, int NCf, int NCor, int NTherm, double a, double L, double seed, Correlator *correlator, Action *S);
    ~Metropolis();
    void latticeSetup(SU3MatrixGenerator *SU3Generator);
    void runMetropolis();
    void getStatistics();
    // Data outputters
//    void writeStatisticsToFile(const char *filename);
    void writeDataToFile(const char *filename);
    void writeConfigurationToFile(std::__1::string filename);
    void loadFieldConfiguration(std::__1::string filename);

    // Setters
    void setAction(Action *S) { m_S = S; }
    void setCorrelator(Correlator *correlator) { m_correlator = correlator; }
    void setN(int N) { m_N = N; }
    void setNT(int N_T) { m_N_T = N_T; }
    void setNCf(int NCf) { m_NCf = NCf; }
    void setEpsilon(double epsilon) { m_epsilon = epsilon; }
    void setUpdateFrequency(int nUpdates) { m_nUpdates = nUpdates; }

    // Getters
    int getN() { return m_N; }
    int getNT() { return m_N; }
    int getNCf() { return m_NCf; }
    int getEpsilon() { return m_epsilon; }

    // Printers
    void printEnergies();
    void printAcceptanceRate();
};

#endif // METROPOLIS_H
