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

System::System()
{
    /*
     * Class for calculating correlators using the System algorithm.
     * Reads in from Parameters, and sets the action accordingly.
     */
    if (Parallel::ParallelParameters::active) {
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
    Parallel::Communicator::setBarrier();
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
    bool topct = false;
    bool energyTopcFieldDensity = false;
    for (unsigned int i = 0; i < obsList.size(); i++) {
        if (obsList[i] == "plaq") plaq = true;
        if (obsList[i] == "topc") topc = true;
        if (obsList[i] == "energy") energy = true;
        if (obsList[i] == "topct") topct = true;
        if (obsList[i] == "energyTopcFieldDensity") energyTopcFieldDensity = true;
    }
    if ((topc || energy) && !topct) {
        // Initializes the full mechinery except for the QxyzQt sampler
        if (flow) {
            m_flowCorrelator = new MasterSampler(flow);
        } else {
            m_correlator = new MasterSampler(flow);
        }
    } else if (plaq && !topc && !energy && !topct) {
        // Initialize plaquette sampler when no other samplers are specified
        if (flow) {
            m_flowCorrelator = new Plaquette(flow);
        } else {
            m_correlator = new Plaquette(flow);
        }
    } else if (energyTopcFieldDensity && !plaq && !topc && !energy && !topct) {
        // Initates energyTopcFieldDensity if only that and no other samplers are specified.
        // Also then sets number of configurations to generate to one,
        // in order to avoid writing out too many fields.
        m_NCf = 1;
        Parameters::setNCf(m_NCf);
        m_writeConfigsToFile = false;
        Parameters::setStoreConfigurations(m_writeConfigsToFile);
        if (flow) {
            m_flowCorrelator = new LatticeActionChargeDensity(flow);
        } else {
            m_correlator = new LatticeActionChargeDensity(flow);
        }
    } else {
        if (flow) {
            m_flowCorrelator = new MasterSamplerTopcXYZ(flow);
        } else {
            m_correlator = new MasterSamplerTopcXYZ(flow);
        }
    }
}


System::~System()
{
    /*
     * Class destructor
     */
    delete m_S;
    delete m_SU3Generator;
    delete [] m_lattice;
    delete m_correlator;

    // Deleting flow related variables - should strictly speaking not be necesseary.
    delete [] m_flowLattice;
    delete m_flowCorrelator;
    delete m_flow;
}

void System::subLatticeSetup()
{
    /*
     * Sets up the sub-lattices.
     */
    Parallel::Communicator::initializeSubLattice();
    m_N = Parameters::getN();
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
        m_flow = new Flow(m_S); // Ensure it does not go out of scope
        m_flowLattice = new Lattice<SU3>[4];
        for (int mu = 0; mu < 4; mu++) {
            m_flowLattice[mu].allocate(m_N);
        }
    }
    IO::FieldIO::init();
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
    if (Parallel::ParallelParameters::active) {
        subLatticeSetup();
        if (!Parameters::getLoadFieldConfigurations() && !Parameters::getLoadConfigAndRun()) {
            if (Parameters::getHotStart()) {
                // All starts with a completely random matrix.
                for (int mu = 0; mu < 4; mu++)
                {
                    for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
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
                    for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
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
    Parallel::Communicator::setBarrier();
}

void System::run()
{
    /*
     * Overarching run function, which sets us of to three posibilities:
     *  - run regular metropolis to generate configurations
     *  - load a thermalized configuration and run metropolis to generate configurations
     *  - loads a set of configurations and flows them
     */
    if (Parallel::ParallelParameters::active) {
        if (!Parameters::getLoadFieldConfigurations() && !Parameters::getLoadConfigAndRun()) {
            // Run regular metropolis
            runMetropolis();
        } else if (Parameters::getLoadConfigAndRun()) {
            // Loads a configuration(which is assumed to be thermalized), and then continue generating from that configuration
            loadConfigurationAndRunMetropolis();
        } else {
            // Run flow on configurations;
            flowConfigurations();
        }
    }
    Parallel::Communicator::setBarrier(); // For waiting for all processors to finish at equal times
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
            printf("\n%-4d ",0);
        }
        m_correlator->printObservable(0);
    }
    // Running thermalization
    for (unsigned int iTherm = 1; iTherm < m_NTherm + 1; iTherm++)
    {
        // Pre update time
        m_preUpdate = steady_clock::now();

        // Pre timer
        update();

        // Post timer
        m_updateTime = duration_cast<duration<double>>(steady_clock::now() - m_preUpdate);
        m_updateStorerTherm += m_updateTime.count();
        if (m_processRank == 0 && iTherm % 20 == 0) { // Progress and Avg. time per update every 10th update
            printf("\n%6.2f %% done. Avg. update time: %10.6f sec", iTherm/double(m_NTherm)*100, m_updateStorerTherm/double(iTherm));
        }

        // Print correlator every somehting or store them all(useful when doing the thermalization).
        if (m_storeThermalizationObservables) {
            // Calculating the correlator
            m_correlator->calculate(m_lattice,iTherm);
            if (m_processRank == 0) printf("\n%-4d ",iTherm);
            m_correlator->printObservable(iTherm);
        }
    }

    // Printing out the avg. update time one more time at the end, to avoid unfinnished percentage sign
    if (m_processRank == 0) {
        printf("\n%6.2f %% done. Avg. update time: %10.6f sec", 100.0, m_updateStorerTherm/double(m_NTherm + 1));
    }

    // Taking the average of the acceptance rate across the processors.
    if (m_NTherm != 0) {
        MPI_Allreduce(&m_acceptanceScore,&m_acceptanceScore,1,MPI_DOUBLE,MPI_SUM,Parallel::ParallelParameters::ACTIVE_COMM);
        // Printing post-thermalization correlator and acceptance rate
        if (m_processRank == 0) {
            printf("\nTermalization complete. Acceptance rate: %f",m_acceptanceScore/double(m_NTherm));
        }
    }
}

void System::updateLink(unsigned int iSite, int mu)
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
                        for (unsigned int n = 0; n < m_NUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
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
    // Updates acceptance storer
    m_acceptanceScore += double(m_acceptanceCounter)/double(4*m_NUpdates*Parameters::getSubLatticeSize());
    m_acceptanceCounter = 0;
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
        m_correlator->printHeader();
        printf(" Avg.Update-time  Accept/reject Config-filename");
    }

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;
    m_acceptanceScore = 0;

    // Variables for checking performance of the update.
    m_updateStorer = 0;

    // Main part of algorithm
    for (unsigned int iConfig = 0; iConfig < m_NCf; iConfig++)
    {
        for (unsigned int i = 0; i < m_NCor; i++) // Updating NCor times before updating the lattice
        {
            // Pre timer
            m_preUpdate = steady_clock::now();
            update();

            // Post timer
            m_updateTime = duration_cast<duration<double>>(steady_clock::now() - m_preUpdate);
            m_updateStorer += m_updateTime.count();
        }

        if (Parameters::getDebug()) {
            Parallel::Communicator::checkLattice(m_lattice, "Configuration is corrupt in system right before copying to flow lattice");
        }

        // Flowing configuration
        if (m_NFlows != 0) {
            // Copys lattice into flow-lattice in order to avoid overwriting when generating configurations
            copyToFlowLattice();

            if (Parameters::getDebug()) {
                Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt in system right after copying to flow lattice");
            }

            flowConfiguration(iConfig);
        }

        // Averaging the observable values. Avoids calculating twice if we are flowing
        if (m_NFlows == 0) {
            m_correlator->calculate(m_lattice,iConfig + m_NThermSteps);
        } else {
            m_correlator->copyObservable(iConfig + m_NThermSteps, m_flowCorrelator->getObservablesVector(0));
        }

        if (m_processRank == 0) {
            // Printing the observables
            printf("\n%-4d ", iConfig);
        }
        m_correlator->printObservable(iConfig + m_NThermSteps);

        if (m_processRank == 0) {
            // Adds time per cfg to output
            printf(" %-12.8f",m_updateStorer/double((iConfig+1)*m_NCor));

            // Adding the acceptance ratio
            if (iConfig % 10 == 0) {
                printf("     %-10.8f", m_acceptanceScore/double((iConfig+1)*m_NCor));
            } else {
                printf("               ");
            }
        }

        if (Parameters::getDebug()) {
            Parallel::Communicator::checkLattice(m_lattice, "Configuration is corrupt in system before write field to file");
        }

        // Writing field config to file
        if (m_writeConfigsToFile)
        {
            IO::FieldIO::writeFieldToFile(m_lattice,iConfig);
        }

    }

    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceScore,&m_acceptanceScore,1,MPI_DOUBLE,MPI_SUM,Parallel::ParallelParameters::ACTIVE_COMM);
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
    m_correlator->writeObservableToFile(getAcceptanceRate()); // Runs statistics, writes to file, and prints results (if verbose is on)
    m_correlator->printStatistics();
}

void System::flowConfiguration(unsigned int iConfig)
{
    /*
     * Flows configuration, performs flow statistics and writes it to a file.
     */

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 0.");
    }

    // After each configuration has been flowed, the values must be resetted.
    m_flowCorrelator->reset();

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 1.");
    }

    // Calculates the flow observables at zeroth flow time
    m_flowCorrelator->calculate(m_flowLattice,0);
    if (Parameters::getVerbose()) {
        m_flowCorrelator->printObservable(0);
    }

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 2.");
    }

    // Runs the flow
    for (unsigned int iFlow = 0; iFlow < m_NFlows; iFlow++)
    {
        if (Parameters::getDebug()) {
            Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 2.5, before flow.");
        }

        m_flow->flowField(m_flowLattice);

        if (Parameters::getDebug()) {
            Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 2.5, after flow.");
        }

        m_flowCorrelator->calculate(m_flowLattice,iFlow + 1);


        if (Parameters::getDebug()) {
            Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 2.5, after flow and after correlator calculation.");
        }

        if (Parameters::getVerbose()) {
            m_flowCorrelator->printObservable(iFlow + 1);
        }
    }


    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 3.");
    }

    // Write flow data to file
    m_flowCorrelator->writeFlowObservablesToFile(iConfig);

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(m_flowLattice, "Configuration is corrupt at flowConfiguration pt 4.");
    }
}

void System::loadConfigurationAndRunMetropolis()
{
    /*
     * Method for loading and running Metropolis updates from a specific configuration.
     */
    m_systemIsThermalized = true;

    // Loads the configuration we are flowing from(should only be one).
    std::vector<std::string> configurationNames = Parameters::getFieldConfigurationFileNames();
    if (configurationNames.size() != 1) {
        Parallel::Communicator::MPIExit("Should only be one configuration loaded(there are currently " + std::to_string(configurationNames.size()) + "configurations.");
    }

    // Loads configuration, either in chroma format(reversed doubles) or regular format.
    if (!Parameters::getLoadChromaConfigurations()) {
        load(configurationNames[0]);
    } else {
        loadChroma(configurationNames[0]);
    }
    runMetropolis();
}


void System::load(std::string configurationName)
{
    /*
     * Method for loading regular a configuration and continuing and evolving it(without the need for any thermalization)
     */
    m_systemIsThermalized = true;
    m_storeThermalizationObservables = false;
    if (m_NFlows != 0 && !Parameters::getLoadConfigAndRun()) {
        if (Parallel::Communicator::getProcessRank() == 0) cout << "\nLoading flow lattice @ line 606" << endl;
        IO::FieldIO::loadFieldConfiguration(configurationName,m_flowLattice);
    } else {
        if (Parallel::Communicator::getProcessRank() == 0) cout << "\nLoading non-flow lattice @ line 609" << endl;
        IO::FieldIO::loadFieldConfiguration(configurationName,m_lattice);
    }
}

void System::loadChroma(std::string configurationName)
{
    /*
     * Method for loading regular a configuration and continuing and evolving it(without the need for any thermalization)
     */
    m_systemIsThermalized = true;
    m_storeThermalizationObservables = false;
    if (m_NFlows != 0 && !Parameters::getLoadConfigAndRun()) {
        IO::FieldIO::loadChromaFieldConfiguration(configurationName,m_flowLattice);
    } else {
        IO::FieldIO::loadChromaFieldConfiguration(configurationName,m_lattice);
    }
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
    return m_acceptanceScore/double(m_NCf*m_NCor*Parallel::Communicator::getNumProc());
}
