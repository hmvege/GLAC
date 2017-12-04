#include "system.h"
#include "config/parameters.h"
#include "parallelization/parallel.h"
#include "io/fieldio.h"
#include <cmath>    // For exp()
#include <mpi.h>

using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

//System::System(Correlator *correlator, Action *S, Flow *F, Correlator *flowCorrelator)
System::System()
{
    /*
     * Class for calculating correlators using the System algorithm.
     * Reads in from Parameters, and sets the action accordingly.
     */
    // Retrieving communication related variables
    m_processRank                       = Parallel::Communicator::getProcessRank();
    // Retrieving program parameters
    m_latticeSize                       = Parameters::getLatticeSize();
    m_NCf                               = Parameters::getNCf();
    m_NCor                              = Parameters::getNCor();
    m_NTherm                            = Parameters::getNTherm();
    m_NUpdates                          = Parameters::getNUpdates();
    m_NFlows                            = Parameters::getNFlows();
    m_storeThermalizationObservables    = Parameters::getStoreThermalizationObservables();
    m_writeConfigsToFile                = Parameters::getStoreConfigurations();
    // Sets pointers to use
    m_SU3Generator                      = new SU3MatrixGenerator;
    // Initializing the Mersenne-Twister19937 RNG for the Metropolis algorithm
    m_generator                         = std::mt19937_64(Parameters::getMetropolisSeed());
    m_uniform_distribution              = std::uniform_real_distribution<double>(0,1);
}

void System::setAction()
{
    /*
     * Only one action available at the moment
     */
    m_S = new WilsonGaugeAction;
}


void System::setObservable(std::vector<std::string> obsList, bool flow)
{
    /*
     * Sets the observables to be sampled.
     * Arguments:
     *  obsList     : vector of string
     *  flow        : flag if we are initializing a flow variable
     */
    bool plaq = false;
    bool topc = false;
    bool energy = false;
    for (unsigned int i = 0; i < obsList.size(); i++) {
        if (obsList[i] == "plaq") plaq = true;
        if (obsList[i] == "topc") topc = true;
        if (obsList[i] == "energy") energy = true;
        // ADD USER OBS?
    }
    if (plaq && !topc && !energy) {
        // Initialize plaquette sampler
        if (flow) {
            m_flowCorrelator = new Plaquette(flow);
        } else {
            m_correlator = new Plaquette(flow);
        }
    }
    else {
        // All other cases, might as well initialize the full machinery
        if (flow) {
            m_flowCorrelator = new ObservableSampler(flow);
        } else {
            m_correlator = new ObservableSampler(flow);
        }
    }
}


System::~System()
{
    /*
     * Class destructor
     */
    delete [] m_lattice;
    if (m_NFlows != 0) delete [] m_flowLattice;
}

void System::subLatticeSetup()
{
    /*
     * Sets up the sub-lattices.
     */
    Parallel::Communicator::initializeSubLattice();
    Parameters::getN(m_N);
    m_subLatticeSize = Parameters::getSubLatticeSize();
    // Creates/allocates (sub) lattice
    m_lattice = new Lattice<SU3>[4];
    for (int mu = 0; mu < 4; mu++) {
        m_lattice[mu].allocate(m_N);
    }
    // Sets pointers
    setAction();
    setObservable(Parameters::getObservablesList(),false);
    // Passes the index handler and dimensionality to the action and correlator classes.
    if (m_NFlows != 0) {
        setObservable(Parameters::getFlowObservablesList(),true);
        m_flowCorrelator->setN(m_N);
        m_flowCorrelator->setLatticeSize(m_subLatticeSize);
        m_flow = new Flow(m_S); // Ensure it does not go out of scope
        m_flowLattice = new Lattice<SU3>[4];
    }
    IO::FieldIO::init();
    m_S->setN(m_N);
    m_correlator->setN(m_N);
    m_correlator->setLatticeSize(m_subLatticeSize);
}

void System::copyToFlowLattice()
{
    /*
     * Small function for copying the lattice to the flow lattice,
     * as we need to ensure the old lattice remains unchanged.
     */
    for (int mu = 0; mu < 4; mu++) {
        m_flowLattice[mu] = m_lattice[mu];
    }
}

void System::latticeSetup()
{
    /*
     * Sets up the lattice and its matrices.
     */
    subLatticeSetup();
    if (!Parameters::getLoadFieldConfigurations()) {
        if (Parameters::getHotStart()) {
            // All starts with a completely random matrix.
            for (int mu = 0; mu < 4; mu++)
            {
                for (int iSite = 0; iSite < m_subLatticeSize; iSite++)
                {
                    if (Parameters::getRSTHotStart())
                    {
                        m_lattice[mu][iSite] = m_SU3Generator->generateRST(); // Random close to unity
                    } else {
                        m_lattice[mu][iSite] = m_SU3Generator->generateRandom(); // Fully random
                    }
                }
            }
        } else {
            // Cold start: everything starts out at unity.
            for (int mu = 0; mu < 4; mu++)
            {
                for (int iSite = 0; iSite < m_subLatticeSize; iSite++)
                {
                    m_lattice[mu][iSite].identity();
                }
            }
        }
        if (m_processRank == 0) printf("\nLattice setup complete\n");
    }
    // Prints system info after setup is complete
    SysPrint::printSystemInfo();
}

void System::run()
{
    if (!Parameters::getLoadFieldConfigurations()) {
        // Run regular metropolis
        runMetropolis();
    } else {
        // Run flow on configurations;
        flowConfigurations();
    }
}

void System::thermalize()
{
    /*
     * Function for thermalizing the system.
     */
    if (m_processRank == 0) printf("\nInitiating thermalization.");
    if (m_storeThermalizationObservables) {
        // Storing the number of shifts that are needed in the observable storage container. Will be used for main for loop shifting.
        m_NThermSteps = 1 + m_NTherm;
        // Calculating correlator before any updates have began.
        m_correlator->calculate(m_lattice,0); // Averaging is done if specified for observable in statistics after run is done.
        if (m_processRank == 0) {
            printf("\ni    ");
            m_correlator->printHeader();
            m_correlator->printObservable(0);
        }
    }
    // Running thermalization
    for (int iTherm = 1; iTherm < m_NTherm + 1; iTherm++)
    {
        // Pre update time
        m_preUpdate = steady_clock::now();

        // Pre timer
        update();

        // Post timer
        m_updateTime = duration_cast<duration<double>>(steady_clock::now() - m_preUpdate);
        m_updateStorerTherm += m_updateTime.count();
        if (m_processRank == 0 && iTherm % 20 == 0) { // Progress and Avg. time per update every 10th update
            printf("\r%6.2f %% done. Avg. update time: %10.6f sec", iTherm/double(m_NTherm)*100, m_updateStorerTherm/double(iTherm));
            if (!m_storeThermalizationObservables) {
                // Flushes if we are not storing any thermalization observables
                std::fflush(stdout);
            }
        }

        // Print correlator every somehting or store them all(useful when doing the thermalization).
        if (m_storeThermalizationObservables) {
            // Calculating the correlator
            m_correlator->calculate(m_lattice,iTherm);
            if (m_processRank == 0) {
                m_correlator->printObservable(iTherm);
            }
        }
    }

    // Printing out the avg. update time one more time at the end, to avoid unfinnished percentage sign
    if (m_processRank == 0) {
        printf("\r%6.2f %% done. Avg. update time: %10.6f sec", 100.0, m_updateStorerTherm/double(m_NTherm + 1));
    }

    // Taking the average of the acceptance rate across the processors.
    if (m_NTherm != 0) {
        MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
        // Printing post-thermalization correlator and acceptance rate
        if (m_processRank == 0) {
            printf("\nTermalization complete. Acceptance rate: %f",double(m_acceptanceCounter)/double(4*m_latticeSize*m_NUpdates*m_NTherm));
        }
    }
}

void System::updateLink(int iSite, int mu)
{
    /*
     * Private function used for updating our system. Updates a single gauge link.
     * Arguments:
     *  i   : spacetime getIndex
     *  mu  : Lorentz getIndex
     */
//    m_updatedMatrix = m_SU3Generator->generateRandom()*m_lattice[latticeIndex].U[mu]; // Shorter method of updating matrix
    m_updatedMatrix = m_SU3Generator->generateRST()*m_lattice[mu][iSite]; // Shorter method of updating matrix
}

void System::update()
{
    /*
     * Sweeps the entire Lattice, and gives every matrix a chance to update.
     */
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice, x, y, z, t, mu);
                        for (int n = 0; n < m_NUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(Parallel::Index::getIndex(x,y,z,t), mu);
                            if (exp(-m_S->getDeltaAction(m_lattice[mu][Parallel::Index::getIndex(x,y,z,t)], m_updatedMatrix)) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[mu][Parallel::Index::getIndex(x,y,z,t)] = m_updatedMatrix;
                                m_acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void System::runMetropolis()
{
    /*
     * Runs the generation of gauge field configurations through the Metropolis algorithm.
     */
    // Variables for checking performance of the thermalization update.
    m_updateStorerTherm = 0;

    // System thermalization
    if (!m_systemIsThermalized) {
        thermalize();
    }

    // Printing header for main run
    if (m_processRank == 0) {
        printf("\ni    ");
        if (m_NFlows == 0) {
            m_correlator->printHeader();
        } else {
            m_flowCorrelator->printHeader();
        }
        printf(" Avg.Update-time  Accept/reject");
    }

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;

    // Variables for checking performance of the update.
    m_updateStorer = 0;

    // Main part of algorithm
    for (int iConfig = 0; iConfig < m_NCf; iConfig++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the lattice
        {
            // Pre timer
            m_preUpdate = steady_clock::now();

            update();

            // Post timer
            m_updateTime = duration_cast<duration<double>>(steady_clock::now() - m_preUpdate);
            m_updateStorer += m_updateTime.count();
        }
        // Flowing configuration
        if (m_NFlows != 0) {
            // Copys lattice into flow-lattice in order to avoid overwriting when generating configurations
            copyToFlowLattice();
            flowConfiguration(iConfig);
        }

        // Averaging the observable values. Avoids calculating twice if we are flowing
        if (m_NFlows == 0) {
            m_correlator->calculate(m_lattice,iConfig + m_NThermSteps);
        } else {
            m_correlator->copyObservable(iConfig, m_flowCorrelator->getObservablesVector(0));
        }

        if (m_processRank == 0) {
            // Printing the observables
            printf("\n%-4d ",iConfig);
            m_correlator->printObservable(iConfig);
            printf(" %-12.8f",m_updateStorer/double((iConfig+1)*m_NCor));
            // Adding the acceptance ratio
            if (iConfig % 10 == 0) {
                printf("     %-12.8f", double(m_acceptanceCounter)/double(4*m_subLatticeSize*(iConfig+1)*m_NUpdates*m_NCor));
            }
        }
        // Writing field config to file
        if (m_writeConfigsToFile)
        {
            IO::FieldIO::writeFieldToFile(m_lattice,iConfig);
        }

    }
    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    if (m_processRank == 0) {
        printf("\n");
        SysPrint::printLine();
        printf("System completed.");
        printf("\nAcceptancerate: %.16f ", getAcceptanceRate());
        printf("\nAverage update time: %.6f sec.", m_updateStorer/double(m_NCf*m_NCor));
        printf("\nTotal update time for %d updates: %.6f sec.\n", m_NCf*m_NCor, m_updateStorer + m_updateStorerTherm);
        SysPrint::printLine();
    }
    m_correlator->runStatistics();
    m_correlator->writeStatisticsToFile(getAcceptanceRate()); // Runs statistics, writes to file, and prints results (if verbose is on)
    m_correlator->printStatistics();
}

void System::flowConfiguration(int iConfig)
{
    /*
     * Flows configuration, performs flow statistics and writes it to a file.
     */
    // After each configuration has been flowed, the values must be resetted.
    m_flowCorrelator->reset();
    // Calculates the flow observables at zeroth flow time
    m_flowCorrelator->calculate(m_flowLattice,0);
    if (Parameters::getVerbose()) {
        m_flowCorrelator->printObservable(0);
    }
    // Runs the flow
    for (int iFlow = 0; iFlow < m_NFlows; iFlow++)
    {
        m_flow->flowField(m_flowLattice);
        m_flowCorrelator->calculate(m_flowLattice,iFlow + 1);
        if (Parameters::getVerbose()) {
            m_flowCorrelator->printObservable(iFlow + 1);
        }
    }
    // Write flow data to file
    m_flowCorrelator->writeFlowObservablesToFile(iConfig);
}


void System::load(std::string configurationName)
{
    /*
     * Method for loading regular a configuration and continuing and evolving it(without the need for any thermalization)
     */
    m_systemIsThermalized = true;
    m_storeThermalizationObservables = false;
    IO::FieldIO::loadFieldConfiguration(configurationName,m_flowLattice);
}

void System::loadChroma(std::string configurationName)
{
    /*
     * Method for loading regular a configuration and continuing and evolving it(without the need for any thermalization)
     */
    m_systemIsThermalized = true;
    m_storeThermalizationObservables = false;
    IO::FieldIO::loadChromaFieldConfiguration(configurationName,m_flowLattice);
}

void System::flowConfigurations()
{
    /*
     * Method for flowing several configurations given as a vector of strings.
     */
    // Loads the vector of configurations to flow.
    std::vector<std::string> configurationNames = Parameters::getFieldConfigurationFileNames();
    for (unsigned int i = 0; i < configurationNames.size(); i++) {
        // Loads configuration, either in chroma format(reversed doubles) or regular format.
        if (!Parameters::getLoadChromaConfigurations()) {
            load(configurationNames[i]);
        } else {
            loadChroma(configurationNames[i]);
        }
        // Prints a new header for each flow.
        if (Parallel::Communicator::getProcessRank() == 0 && Parameters::getVerbose()) {
            m_flowCorrelator->printHeader();
        }
        // Flows the configuration loaded
        flowConfiguration(i);
    }
    if (m_processRank==0) printf("\nFlowing of %lu configurations done.", configurationNames.size());
}

double System::getAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the System algorithm.
     */
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_NUpdates*m_latticeSize*4); // Times 4 from the Lorentz indices
}
