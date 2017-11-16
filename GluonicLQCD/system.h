#ifndef SYSTEM_H
#define SYSTEM_H

#include <random>
#include <chrono>
#include "actions/action.h"
#include "observables/correlator.h"
#include "config/parameters.h"
#include "math/matrices/su3matrixgenerator.h"
#include "math/latticemath.h"

#include "parallelization/neighbourlist.h" // Make this into
#include "parallelization/neighbours.h"
#include "parallelization/index.h"

#include "flow/flow.h"
#include "io/fieldio.h"
#include "io/observablesio.h"

#include "observables/observablesampler.h"
//#include "observablestorer.h"

using std::chrono::steady_clock;
using std::chrono::duration;

class Action;

class System
{
private:
    // Lattice sizes
    int m_NSpatial;
    int m_NTemporal;
    unsigned int m_N[4];

    // Beta value constant
    double m_beta;

    // Variable for storing the thermalization observables
    bool m_storeThermalizationObservables = false;
    bool m_systemIsThermalized = false;

    // Updating constants
    int m_NCf;
    int m_NCor;
    int m_NTherm;
    int m_NUpdates; // N updates before calculating the action, as that is costly
    int m_NFlows;

    // Variable for storing how many steps we are shifting in the observables storage array if we choose to store the thermalization variables
    int m_NThermSteps = 0;

    // For handling the acceptance rate
    unsigned long int m_acceptanceCounter = 0;
    double getAcceptanceRate();

    // Variable for storing if the sub lattice size has been preset.
    bool m_subLatticeSizePreset = false;

    // For have the possibility to start lattice with SU3 RST random, and not fully random SU3
    bool m_RSTInit = false;

    // Paralellization setup
    int m_numprocs;
    int m_processRank;
    int m_processorsPerDimension[4];
    int m_subLatticeSize;
    void subLatticeSetup();
    int linkDoubles = 72;
    int linkSize = linkDoubles*sizeof(double);

    // Lattice variables
    int m_latticeSize;
    Links * m_lattice;
    SU3 m_updatedMatrix;

    // Parallelization variables and functions
    Neighbours * m_neighbourLists = nullptr;

    // Variables used to perform statistics
    double * m_observable;
    double * m_observablePreThermalization;
    double * m_observableSquared;
    double m_averagedObservable = 0; // Change these to not have m_ convention
    double m_varianceObservable = 0;
    double m_stdObservable = 0;

    // Time counting
    steady_clock::time_point m_preUpdate, m_postUpdate;
    duration<double> m_updateTime;
    double m_updateStorer = 0;
    double m_updateStorerTherm = 0;

    // Storing the action as a pointer
    Action *m_S = nullptr;

    // Config correlator
    Correlator * m_correlator = nullptr;

    // Flow correlator
    Correlator * m_flowCorrelator = nullptr;

    // Flow
    Flow * m_Flow = nullptr;
    void flowConfiguration(int iConfig);

    // Function for updating our system using the Metropolis algorithm
    void update();
    void updateLink(int latticeIndex, int mu);

    // Thermalization function
    void thermalize();

    // Input/output locations REDUNDANT
    std::string m_pwd = "";
    std::string m_batchName = "";
    std::string m_inputFolder = "/input/";
    std::string m_outputFolder = "/output/"; // On mac, do not need ../

    // SU3 generator
    SU3MatrixGenerator *m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;

    inline void printLine();
public:
    System(double seed, Correlator *correlator, Action *S, Flow *F, Correlator *flowCorrelator);
    ~System();
    void runMetropolis(bool storeThermalizationObservables, bool writeConfigsToFile);
    void latticeSetup(SU3MatrixGenerator *SU3Generator, bool hotStart);
    void runBasicStatistics();

    // Functions loading fields configurations from file
    void loadChroma(std::string configurationName);
    void load(std::string configurationName);
    void flowConfigurations(std::vector<std::string> configurationNames);

    // Object setters
    void setAction(Action *S) { m_S = S; }
    void setCorrelator(Correlator *correlator) { m_correlator = correlator; }
    void setFlowCorrelator(Correlator *flowCorrelator) { m_flowCorrelator = flowCorrelator; }
    void setFlow(Flow *F) { m_Flow = F; }
    void setSU3ExpFunc(SU3Exp * SU3ExpFunc) { m_Flow->setSU3ExpFunc(SU3ExpFunc); }

    // Variable setters
    void setSubLatticeDimensions(int *NSub);
    void setLatticeInitRST(bool RSTInit) { m_RSTInit = RSTInit; }

    // Printers
    void printRunInfo(bool verbose);
    void printAcceptanceRate();
};

#endif // METROPOLIS_H
