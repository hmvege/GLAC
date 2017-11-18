#include <random>   // For Mersenne-Twister19937
#include <chrono>
#include <cmath>    // For exp()
#include <fstream> // REDUNDANT
#include <iostream> // REDUNDANT
#include <iomanip> // REDUNDANT
#include <cstdio>   // For io C-style handling. // REDUNDANT
#include <cstdlib> // REDUNDANT
#include <mpi.h>
#include "system.h"

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
     * Takes an action object as well as a Gamma functional to be used in the action.
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
    m_SU3Generator = new SU3MatrixGenerator;
    setAction();
    setObservable(Parameters::getObservablesList(),false);
    if (m_NFlows != 0) {
        setObservable(Parameters::getFlowObservablesList(),true);
        m_flow = new Flow(m_S); // Ensure it does not go out of scope
    }

    // Initializing the Mersenne-Twister19937 RNG for the Metropolis algorithm
    m_generator = std::mt19937_64(Parameters::getMetropolisSeed());
    m_uniform_distribution = std::uniform_real_distribution<double>(0,1);
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
    bool plaq = false;
    bool topc = false;
    bool energy = false;
    for (unsigned int i = 0; i < obsList.size(); i++) {
        if (obsList[i] == "plaquette") plaq = true;
        if (obsList[i] == "topc") topc = true;
        if (obsList[i] == "energy") energy = true;
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
}

void System::subLatticeSetup()
{
    /*
     * Sets up the sub-lattices.
     */
    Parallel::Communicator::initializeSubLattice();
    Parameters::getN(m_N);
    m_subLatticeSize = Parameters::getSubLatticeSize();
    // Creates (sub) lattice
    m_lattice = new Links[m_subLatticeSize];
    // Passes the index handler and dimensionality to the action and correlator classes.
    m_S->setN(m_N);
    m_correlator->setN(m_N);
    m_correlator->setLatticeSize(m_subLatticeSize);
}

void System::latticeSetup()
{
    /*
     * Sets up the lattice and its matrices.
     */
    subLatticeSetup();
    if (Parameters::getHotStart()) {
        // All starts with a completely random matrix.
        for (int i = 0; i < m_subLatticeSize; i++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                if (Parameters::getRSTInit())
                {
                    m_lattice[i].U[mu] = m_SU3Generator->generateRST(); // Random close to unity
                } else {
                    m_lattice[i].U[mu] = m_SU3Generator->generateRandom(); // Fully random
                }
            }
        }
    } else {
        // Cold start: everything starts out at unity.
        for (int i = 0; i < m_subLatticeSize; i++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                m_lattice[i].U[mu].identity();
            }
        }
    }
    if (m_processRank == 0) {
        printf("\nLattice setup complete");
    }
}

void System::thermalize()
{
    /*
     * Function for thermalizing the system.
     */
    if (m_processRank == 0) printf("\nInitiating thermalization.");
    if (m_storeThermalizationObservables) {
        // Storing the number of shifts that are needed in the observable storage container.
        m_NThermSteps = 1 + m_NTherm;
        // Calculating correlator before any updates have began.
        m_correlator->calculate(m_lattice,0);
//        // Summing and sharing correlator to all processors before any updates has begun
//        MPI_Allreduce(&m_observablePreThermalization[0], &m_observablePreThermalization[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        // Dividing by the number of processors in order to get the correlator.
//        m_observablePreThermalization[0] /= double(m_numprocs);
        if (m_processRank == 0) {
            printf("\ni    Observable   ");
//            printf("\n%-4d %-12.8f",0,m_observablePreThermalization[0]);
            printf("\n%-4d %-12.8f",0,m_correlator->getObservable(0));
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
        m_postUpdate = steady_clock::now(); // REDUNDANT?
        m_updateTime = duration_cast<duration<double>>(m_postUpdate - m_preUpdate);
        m_updateStorerTherm += m_updateTime.count();
        if (m_processRank == 0) {
            printf("\r%6.2f %% done. ", iTherm/double(m_NTherm));
            if (iTherm % 20 == 0) { // Avg. time per update every 10th update
                printf("Avgerage update time(every 10th): %10.6f sec.", m_updateStorerTherm/double(iTherm));
            }
        }

        // Print correlator every somehting or store them all(useful when doing the thermalization).
        if (m_storeThermalizationObservables) {
            // Calculating the correlator
//            m_observablePreThermalization[i] = m_correlator->calculate(m_lattice);
            m_correlator->calculate(m_lattice,iTherm);

            // Summing and sharing results across the processors
//            MPI_Allreduce(&m_observablePreThermalization[i], &m_observablePreThermalization[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Turn off!
            // Averaging the results
//            m_observablePreThermalization[i] /= double(m_numprocs);

//            if (m_processRank == 0) {
//                printf("\n%-4d %-12.8f",i,m_observablePreThermalization[i]);
//            }
            if (m_processRank == 0) {
                m_correlator->printObservable(iTherm);
//                printf("\n%-4d %-12.8f",iTherm,m_correlator->getObservable(iTherm)); // returns the observable at i(not averaged between by processors?)
            }
        }
    }

    // Taking the average of the acceptance rate across the processors.
    if (m_NTherm != 0) MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    // Printing post-thermalization correlator and acceptance rate
    if (m_processRank == 0 && m_NTherm != 0) printf("\nTermalization complete. Acceptance rate: %f",double(m_acceptanceCounter)/double(4*m_latticeSize*m_NUpdates*m_NTherm));
}

void System::updateLink(int latticeIndex, int mu)
{
    /*
     * Private function used for updating our system. Updates a single gauge link.
     * Arguments:
     *  i   : spacetime getIndex
     *  mu  : Lorentz getIndex
     */
//    m_updatedMatrix = m_SU3Generator->generateRandom()*m_lattice[latticeIndex].U[mu]; // Shorter method of updating matrix
    m_updatedMatrix = m_SU3Generator->generateRST()*m_lattice[latticeIndex].U[mu]; // Shorter method of updating matrix
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
                            if (exp(-m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu)) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[Parallel::Index::getIndex(x,y,z,t)].U[mu] = m_updatedMatrix;
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
    //// TESTS ==============================================================================
//    // Common files
//    loadFieldConfiguration("FlowTestRun_beta6.000000_spatial16_temporal16_threads8_config0.bin"); // 0.59486412, MAC
////    loadFieldConfiguration("msg01.rec02.ildg-binary-data"); // jack
//    double corr = m_correlator->calculate(m_lattice);
//    MPI_Allreduce(&corr, &corr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    corr /= double(m_numprocs);
//    if (m_processRank == 0) cout << "Plaquette value: " << corr << endl;
//    Flow WFlow(m_N, m_beta, m_numprocs, m_processRank);
//    WFlow.setAction(m_S);
//    /// OLD, run tests and compare times?
////    Clover Clov;
////    Clov.initializeIndexHandler(m_indexHandler);
////    Clov.setN(m_N);
////    Clov.setLatticeSize(m_latticeSize);
////    TopologicalCharge TopCharge;
////    TopCharge.initializeIndexHandler(m_indexHandler);
////    TopCharge.setLatticeSize(m_latticeSize);
////    TopCharge.setN(m_N);
////    EnergyDensity Energy(0.0931, m_latticeSize);
////    Energy.initializeIndexHandler(m_indexHandler);

//    ObservableSampler OSampler(m_N,m_subLatticeSize,0.0931);

//    int NFlows = 20;
//    double * m_observableFlow = new double[NFlows];
//    double * m_topologicalCharge = new double[NFlows];
//    double * m_topologicalSusceptibility = new double[NFlows];
//    double * m_actionDensity = new double[NFlows];
//    for (int tau = 0; tau < NFlows; tau++) {
//        m_topologicalCharge[tau] = 0;
//        m_topologicalSusceptibility[tau] = 0;
//        m_observableFlow[tau] = 0;
//        m_actionDensity[tau] = 0;
//    }
//    double updateTime = 0;
//    for (int tau = 0; tau < NFlows; tau++) {
//        m_preUpdate = steady_clock::now();
//        WFlow.flowField(m_lattice);

//        OSampler.calculate(m_lattice);
//        m_observableFlow[tau] = OSampler.getPlaquette();
//        m_topologicalCharge[tau] = OSampler.getTopologicalCharge();
//        m_actionDensity[tau] = OSampler.getEnergyDensity();

//       /// OLD
////        for (unsigned int x = 0; x < m_N[0]; x++) { // CLEAN UP AND MOVE THIS PART INTO ITS OWN CLASS FOR CALCULATING TOP CHARGE AND ENERGY?!
////            for (unsigned int y = 0; y < m_N[1]; y++) { // HIDE IT, AS IT IS BIG AND UGLY!
////                for (unsigned int z = 0; z < m_N[2]; z++) {
////                    for (unsigned int t = 0; t < m_N[3]; t++) {
////                        Clov.calculateClover(m_lattice,x,y,z,t);
////                        m_topologicalCharge[tau] += TopCharge.calculate(Clov.m_clovers);
////                        m_actionDensity[tau] += Energy.calculate(Clov.m_clovers);
////                        m_observableFlow[tau] += m_correlator->calculate(Clov.m_plaquettes);
////                    }
////                }
////            }
////        }

//        MPI_Allreduce(&m_topologicalCharge[tau], &m_topologicalCharge[tau], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&m_actionDensity[tau], &m_actionDensity[tau], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&m_observableFlow[tau], &m_observableFlow[tau], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        m_topologicalSusceptibility[tau] = pow(m_topologicalCharge[tau]*m_topologicalCharge[tau],0.25) * 0.1973/(0.0931*16);
//        m_observableFlow[tau] /= double(m_numprocs);
//        if (m_processRank == 0) printf("\n%5d %-5.4f %-18.16f %-18.16f %-18.16f %-18.16f", tau, 0.0931*sqrt(8*double(0.01*tau)), m_observableFlow[tau], m_topologicalCharge[tau], m_topologicalSusceptibility[tau], m_actionDensity[tau]);

//        updateTime += (duration_cast<duration<double>>(steady_clock::now() - m_preUpdate)).count();

//        if (m_processRank == 0) printf("  Update time: : %-.4f",updateTime / (tau+1));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (m_processRank == 0) printf("\nTime used to flow: %-.4f",updateTime);

//    delete [] m_observableFlow;
//    delete [] m_topologicalCharge;
//    delete [] m_actionDensity;
//    MPI_Finalize(); exit(1);
    //// ===================================================================================
    // Variables for checking performance of the thermalization update.
    m_updateStorerTherm = 0;

    // System thermalization
    if (!m_systemIsThermalized) {
        thermalize();
    }

    // Printing header for main run
    if (m_processRank == 0) {
        printf("\ni    %-*s Avg.Update-time  Accept/reject", m_correlator->getHeaderWidth(),m_correlator->getObservableName().c_str());
    }

    // Setting the System acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;

    // Variables for checking performance of the update.
    m_updateStorer = 0;

    // Main part of algorithm
    for (int iConfig = 0; iConfig < m_NCf; iConfig++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            // Pre timer
            m_preUpdate = steady_clock::now();

            update();

            // Post timer
            m_postUpdate = steady_clock::now();
            m_updateTime = duration_cast<duration<double>>(m_postUpdate - m_preUpdate);
            m_updateStorer += m_updateTime.count();
        }
        // Flowing configuration
        if (m_NFlows != 0) flowConfiguration(iConfig);

        // Averaging the gamma values
        m_correlator->calculate(m_lattice,iConfig);
//        m_observable[iConfig] = m_correlator->calculate(m_lattice);
//        MPI_Allreduce(&m_observable[iConfig], &m_observable[iConfig], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        m_observable[iConfig] /= double(m_numprocs);

        if (m_processRank == 0) {
            // Printing plaquette value
            m_correlator->printObservable(iConfig);
            printf(" %-12.8f",m_updateStorer/double((iConfig+1)*m_NCor));
            // Adding the acceptance ratio
            if (iConfig % 10 == 0) {
                printf("     %-12.8f", double(m_acceptanceCounter)/double(4*m_subLatticeSize*(iConfig+1)*m_NUpdates*m_NCor));
            }
        }
        // Writing field config to file
        if (m_writeConfigsToFile) IO::FieldIO::writeFieldToFile(m_lattice,iConfig);
    }
    // Taking the average of the acceptance rate across the processors.
    MPI_Allreduce(&m_acceptanceCounter,&m_acceptanceCounter,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    if (m_processRank == 0) {
        SysPrint::printLine();
        printf("System completed.");
        printf("\nAcceptancerate: %.16f ", getAcceptanceRate());
        printf("\nAverage update time: %.6f sec.", m_updateStorer/double(m_NCf*m_NCor));
        printf("\nTotal update time for %d updates: %.6f sec.\n", m_NCf*m_NCor, m_updateStorer + m_updateStorerTherm);
        SysPrint::printLine();
    }
    m_correlator->runStatistics();
    m_correlator->writeStatisticsToFile(); // Runs statistics, writes to file, and prints results (if verbose is on)
}

void System::flowConfiguration(int iConfig)
{
    /*
     * Flows configuration, performs flow statistics and writes it to a file.
     */
    for (int iFlow = 0; iFlow < m_NFlows; iFlow++)
    {
        m_flow->flowField(m_lattice);
        m_flowCorrelator->calculate(m_lattice,iFlow + m_NThermSteps);
    }
    /* Make flow statistics, that is, the only stats that is needed is the sum of all configurations, which cant be reached untill end.
         * So, either print flow during the run, or nothing.
         * Then, write the flow values to file. */
    m_flowCorrelator->runStatistics();
    // Write flow data to file
    m_flowCorrelator->writeStatisticsToFile(iConfig);
}


void System::load(std::string configurationName)
{
    /*
     * Method for loading regular a configuration and continuing and evolving it(without the need for any thermalization)
     */
    m_systemIsThermalized = true;
    m_storeThermalizationObservables = false;
    IO::FieldIO::loadFieldConfiguration(configurationName,m_lattice);
}

void System::loadChroma(std::string configurationName)
{
    /*
     * Method for loading regular a configuration and continuing and evolving it(without the need for any thermalization)
     */
    m_systemIsThermalized = true;
    m_storeThermalizationObservables = false;
//    if (Parameters.getConfigType == "chroma") SWITCH TO THIS!
    IO::FieldIO::loadChromaFieldConfiguration(configurationName,m_lattice);
}

void System::flowConfigurations(std::vector<std::string> configurationNames)
{
    /*
     * Method for flowing several configurations given as a vector of strings.
     */
    for (unsigned int i = 0; i < configurationNames.size(); i++) {
        load(configurationNames[i]);
        flowConfiguration(i);
    }
    printf("\nFlowing of %lu configurations done.", configurationNames.size());
}

//void System::runBasicStatistics()
//{
//    /*
//     * Class instance for sampling statistics from our system.
//     */
//    double averagedObservableSquared = 0;
//    // Performing an average over the Monte Carlo obtained values
//    for (int alpha = 0; alpha < m_NCf; alpha++)
//    {
//        m_averagedObservable += m_observable[alpha];
//        averagedObservableSquared += m_observable[alpha]*m_observable[alpha];
//    }
//    m_averagedObservable /= double(m_NCf);
//    averagedObservableSquared /= double(m_NCf);
//    m_varianceObservable = (averagedObservableSquared - m_averagedObservable*m_averagedObservable)/double(m_NCf);
//    m_stdObservable = sqrt(m_varianceObservable);
//    if (m_processRank == 0) {
//        SysPrint::printLine();
//        cout << "Average plaqutte:      " << m_averagedObservable << endl;
//        cout << "Standard deviation:    " << m_stdObservable << endl;
//        cout << "Variance:              " << m_varianceObservable << endl;
//        SysPrint::printLine();
//    }
//}

double System::getAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the System algorithm.
     */
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_NUpdates*m_latticeSize*4); // Times 4 from the Lorentz indices
}
