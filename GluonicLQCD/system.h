#ifndef SYSTEM_H
#define SYSTEM_H

#include <random>
#include "actions/action.h"
#include "correlators/correlator.h"
#include "links.h"
#include "matrices/su3matrixgenerator.h"
#include "parallelization/neighbours.h"
#include "parallelization/indexorganiser.h"

class Action;

class System
{
private:
    // Lattice sizes
    int m_NSpatial;
    int m_NTemporal;
    int m_NTrue[4];
    int m_N[4]; // With phases
    // Updating constants
    int m_NCf;
    int m_NCor;
    int m_NTherm;
    int m_nUpdates = 10; // N updates before calculating the action, as that is costly
    double m_epsilon;
    // For handling the acceptance rate
    int m_acceptanceCounter = 0;
    double getAcceptanceRate();

    // Paralellization setup
    int m_numprocs;
    int m_processRank;
    int m_processorsPerDimension[4];
    void subLatticeSetup();
    IndexOrganiser * m_indexHandler;
    void setIndexHandler(IndexOrganiser *indexHandler);

    // Lattice variables
    int m_latticeSize;
    int m_subLatticeSize;
    int m_trueSubLatticeSize; // With phases
    Links * m_lattice;
    SU3 m_updatedMatrix;

    // Parallelization variables and functions
    Neighbours * m_neighbourLists = nullptr;
    void share();
    void shareFaces();
    void shareEdges();

    // Variables used to perform statistics
    double * m_Gamma;
    double * m_GammaPreThermalization;
    double * m_GammaSquared;
    double m_averagedGamma = 0; // Change these to not have m_ convention
    double m_varianceGamma = 0;
    double m_stdGamma = 0;

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

    // Input/output locations
    std::string m_inputFolder = "../input/";
    std::string m_outputFolder = "../output/";

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;
public:
    System(int NSpatial, int NTemporal, int NCf, int NCor, int NTherm, double seed, Correlator *correlator, Action *S, int numprocs, int processRank);
    ~System();
    void runMetropolis(bool storePreObservables);
    void latticeSetup(SU3MatrixGenerator *SU3Generator, bool hotStart);
    void runBasicStatistics();

    // Data outputters
    void writeDataToFile(std::string filename, bool preThermalizationGamma = true);
    void writeConfigurationToFile(std::string filename);
    void loadFieldConfiguration(std::string filename);

    // Setters
    void setAction(Action *S) { m_S = S; }
    void setCorrelator(Correlator *correlator) { m_correlator = correlator; }
    void setN(int NSpatial) { m_NSpatial = NSpatial; }
    void setNT(int NTemporal) { m_NTemporal = NTemporal; }
    void setNCf(int NCf) { m_NCf = NCf; }
    void setEpsilon(double epsilon) { m_epsilon = epsilon; }
    void setUpdateFrequency(int nUpdates) { m_nUpdates = nUpdates; }

    // Getters
    int getNSpatial() { return m_NSpatial; }
    int getNNTemporal() { return m_NTemporal; }
    int getNCf() { return m_NCf; }
    int getEpsilon() { return m_epsilon; }

    // Printers
    void printEnergies();
    void printAcceptanceRate();
};

#endif // METROPOLIS_H
