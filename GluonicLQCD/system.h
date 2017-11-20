#ifndef SYSTEM_H
#define SYSTEM_H

#include <random>
#include <chrono>
#include "actions/action.h"
#include "actions/wilsongaugeaction.h"
#include "observables/correlator.h"
#include "config/parameters.h"
#include "math/matrices/su3matrixgenerator.h"
#include "math/latticemath.h"

#include "parallelization/neighbourlist.h" // Make this into header?
#include "parallelization/neighbours.h"
#include "parallelization/index.h"

#include "flow/flow.h"
#include "io/fieldio.h"
#include "io/observablesio.h"

#include "observables/observablesampler.h"

using std::chrono::steady_clock;
using std::chrono::duration;

class Action;

class System
{
private:
    /////////////////////////////////////////
    /// Parameters that must be retrieved ///
    /////////////////////////////////////////
    // Lattice sizes
    unsigned int m_N[4];
    int m_latticeSize;
    // Updating constants
    int m_NCf;
    int m_NCor;
    int m_NTherm;
    int m_NUpdates; // N updates before calculating the action, as that is costly
    int m_NFlows;
    // Variable for storing the thermalization observables
    bool m_storeThermalizationObservables = false;
    bool m_systemIsThermalized = false;
    bool m_writeConfigsToFile = false;
    // Paralellization setup
    int m_processRank;
    int m_subLatticeSize;

    /////////////////////////////////////////
    /////// Parameters and functions ////////
    /////////////////////////////////////////
    // Variable for storing how many steps we are shifting in the observables storage array if we choose to store the thermalization variables
    int m_NThermSteps = 0;

    // For handling the acceptance rate
    unsigned long int m_acceptanceCounter = 0;
    double getAcceptanceRate();

    // Function for initializing the sub lattice.
    void subLatticeSetup();

    // Lattice variables
    Links * m_lattice;
    SU3 m_updatedMatrix;

    // Time counting
    steady_clock::time_point m_preUpdate, m_postUpdate;
    duration<double> m_updateTime;
    double m_updateStorer = 0;
    double m_updateStorerTherm = 0;

    // Function for choosing and setting the correlators/observables
    void setObservable(std::vector<std::string> obsList, bool flow);

    // Storing the action as a pointer
    void setAction();
    Action * m_S = nullptr;

    // Config correlator
    Correlator * m_correlator = nullptr;

    // Flow correlator
    Correlator * m_flowCorrelator = nullptr;

    // Flow
    Flow * m_flow = nullptr;
    void flowConfiguration(int iConfig);
    void copyToFlowLattice();
    Links * m_flowLattice;

    // Function for updating our system using the Metropolis algorithm
    void update();
    void updateLink(int latticeIndex, int mu);

    // Thermalization function
    void thermalize();

    // SU3 generator
    SU3MatrixGenerator * m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;
public:
    System();
    ~System();
    void runMetropolis();
    void latticeSetup();

    // Functions loading fields configurations from file
    void loadChroma(std::string configurationName);
    void load(std::string configurationName);
    void flowConfigurations(std::vector<std::string> configurationNames);
};

#endif // METROPOLIS_H
