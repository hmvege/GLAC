/*!
 * \class System
 *
 * \brief System is the class that ties the program together. It initiates and sets up the sub-lattices, and performs holds the Metropolis algorithm.
 *
 * System takes on parameters as they are given by the class Parameters.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef SYSTEM_H
#define SYSTEM_H

#include <random>
#include <chrono>
#include "actions/actions.h"
#include "observables/observables.h"
#include "math/matrices/su3matrixgenerator.h"
#include "math/lattice.h"
#include "flow/flow.h"

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
    std::vector<unsigned int> m_N;
    unsigned int m_latticeSize;
    // Updating constants
    unsigned int m_NCf;
    unsigned int m_NCor;
    unsigned int m_NTherm;
    unsigned int m_NUpdates; // N updates before calculating the action, as that is costly
    unsigned int m_NFlows;
    // Variable for storing the thermalization observables
    bool m_storeThermalizationObservables = false;
    bool m_systemIsThermalized = false;
    bool m_writeConfigsToFile = false;
    // Paralellization setup
    int m_processRank; // Move to communicator/printer...?
    unsigned int m_subLatticeSize;

    /////////////////////////////////////////
    /////// Parameters and functions ////////
    /////////////////////////////////////////
    // Variable for storing how many steps we are shifting in the observables storage array if we choose to store the thermalization variables
    unsigned int m_NThermSteps = 0;

    // For handling the acceptance rate
    unsigned long long int m_acceptanceCounter = 0;
    double m_acceptanceScore = 0;
    double getAcceptanceRate();

    // Function for initializing the sub lattice.
    void subLatticeSetup();

    // Lattice variables
    Lattice<SU3> * m_lattice;
    SU3 m_updatedMatrix;

    // Time counting
    steady_clock::time_point m_preUpdate;
    duration<double> m_updateTime;
    double m_updateStorer = 0;
    double m_updateStorerTherm = 0;

    // Function for choosing and setting the correlators/observables
    void setObservable(const std::vector<std::string> &obsList, const bool flow);

    // Storing the action as a pointer
    void setAction();
    Action * m_S = nullptr;

    // Config correlator
    Correlator * m_correlator = nullptr;

    // Flow correlator
    Correlator * m_flowCorrelator = nullptr;

    // Flow
    Flow * m_flow = nullptr;
    void flowConfiguration(const unsigned int iConfig);
    void copyToFlowLattice();
    Lattice<SU3> * m_flowLattice;

    // Function for updating our system using the Metropolis algorithm
    void update();
    inline void updateLink(const unsigned int iSite, const int mu);

    // Thermalization function
    void thermalize();

    // SU3 generator
    SU3MatrixGenerator * m_SU3Generator = nullptr;

    // RNGs
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_uniform_distribution;

    // Functions loading fields configurations from file
    void loadChroma(const std::string &configurationName);
    void load(const std::string &configurationName);
    void flowConfigurations();
    void loadConfigurationAndRunMetropolis();

    // Function for running metropolis algorithm
    void runMetropolis();
public:
    System();
    ~System();
    void run();
    void latticeSetup();

};

#endif // METROPOLIS_H
