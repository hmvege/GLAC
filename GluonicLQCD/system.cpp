#include <random>   // For Mersenne-Twister19937
#include <chrono>
#include <cmath>    // For exp()
#include <fstream>
#include <iostream>
#include <cstdio>   // For io C-style handling.
#include <cstdlib>
#include "system.h"
#include "actions/action.h"
#include "correlators/correlator.h"
#include "functions.h"
#include "links.h"
#include "matrices/su3matrixgenerator.h"

//TEMP
#include "unittests.h"

using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;


System::System(int N, int N_T, int NCf, int NCor, int NTherm, double seed, Correlator *correlator, Action *S)
{
    /*
     * Class for calculating correlators using the metropolis algorithm.
     * Takes an action object as well as a Gamma functional to be used in the action.
     */
    m_N = N; // Spatial dimensions
    m_N_T = N_T; // Time dimensions
    m_latticeSize = N*N*N*N_T;
    m_NCf = NCf; // Number of configurations to run for
    m_NCor = NCor;
    m_NTherm = NTherm;
    setAction(S);
    setCorrelator(correlator);
    m_lattice = new Links[m_latticeSize]; // Lattice, contigious memory allocation
    m_GammaPreThermalization = new double[m_NTherm+1];
    m_Gamma = new double[m_NCf]; // Correlator values
    m_GammaSquared = new double[m_NCf];

    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uni_dist(0,1);
    m_generator = gen;
    m_uniform_distribution = uni_dist;
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        m_Gamma[alpha] = 0;
        m_GammaSquared[alpha] = 0;
    }
}

System::~System()
{
    /*
     * Class destructor
     */
    delete [] m_lattice;
    delete [] m_Gamma;
    delete [] m_GammaSquared;
}

void System::latticeSetup(SU3MatrixGenerator *SU3Generator)
{
    /*
     * Sets up the lattice and its matrices.
     */
    m_SU3Generator = SU3Generator;
    for (int i = 0; i < m_latticeSize; i++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
//            m_lattice[i].U[mu] = m_SU3Generator->generateRandom();
//            m_lattice[i].U[mu] = m_SU3Generator->generateRST();
            m_lattice[i].U[mu].identity(); // GENERATES IDENTITY FOR TEST! ONE OBSERVABLE SHOULD EQUAL 1!!
        }
    }
}

void System::updateLink(int latticeIndex, int mu)
{
    /*
     * Private function used for updating our system. Updates a single gauge link.
     * Arguments:
     *  i   : spacetime index
     *  mu  : Lorentz index
     */
//    SU3 X = m_SU3Generator->generateRandom(); // Generates a random matrix, SHOULD BE MODIFIED TO X = RST, page 83 Gattinger & Lang
//    SU3 X = m_SU3Generator->generateRST();
//    m_updatedMatrix = X*m_lattice[latticeIndex].U[mu];
    m_updatedMatrix = m_SU3Generator->generateRST()*m_lattice[latticeIndex].U[mu]; // Shorter method of updating matrix
}

void System::update()
{
    /*
     * Sweeps the entire Lattice, and gives every matrix a chance to update.
     */
    for (int x = 0; x < m_N; x++) {
        for (int y = 0; y < m_N; y++) {
            for (int z = 0; z < m_N; z++) {
                for (int t = 0; t < m_N_T; t++) {
                    for (int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice, x, y, z, t, mu);
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(index(x, y, z, t, m_N, m_N_T), mu);
                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, x, y, z, t, mu);
                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[index(x, y, z, t, m_N, m_N_T)].U[mu].copy(m_updatedMatrix);
                                m_acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
}


void System::runMetropolis(bool storeThermalizationObservables)
{
    // TESTS ==============================================================================================
//    loadFieldConfiguration("scalar16cubed16run1");
//    loadFieldConfiguration("parallel8core16cube16");
//    writeConfigurationToFile("unityScalar");
//    loadFieldConfiguration("parallel32core16cube_plaquette0594052");
//    loadFieldConfiguration("unityScalar");
//    m_GammaPreThermalization[0] = m_correlator->calculate(m_lattice);
//    cout << "Pre-thermialization correlator:  " << m_GammaPreThermalization[0] << endl;
//    exit(1);
    // ====================================================================================================
    m_storeThermalizationObservables = storeThermalizationObservables;

    // Timer variables
    steady_clock::time_point preUpdate, postUpdate;
    duration<double> updateTime;
    double updateStorer = 0;

    // Running thermalization
    for (int i = 0; i < m_NTherm; i++)
    {
        preUpdate = steady_clock::now();

        update();

        postUpdate = steady_clock::now();
        updateTime = duration_cast<duration<double>>(postUpdate-preUpdate);
        updateStorer += updateTime.count();

        if (m_storeThermalizationObservables) {
            m_GammaPreThermalization[i+1] = m_correlator->calculate(m_lattice);
            cout << i << " Correlator: " << m_GammaPreThermalization[i+1] << endl;
        }

        if (i % 20 == 0) { // Avg. time per update every 20th update
            cout << "Avgerage update time: " << updateStorer/double(i+1) << " sec" << endl;
        }
    }
    if (m_storeThermalizationObservables) cout << "Post-thermialization correlator: " << m_GammaPreThermalization[m_NTherm] << endl;
    cout << "Termalization complete. Acceptance rate: " << m_acceptanceCounter/double(4*m_latticeSize*m_nUpdates*m_NTherm) << endl;

    // Setting the metropolis acceptance counter to 0 in order not to count the thermalization
    m_acceptanceCounter = 0;

    // Main part of algorithm
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
        }
        m_Gamma[alpha] = m_correlator->calculate(m_lattice);
    }
    cout << "Metropolis completed." << endl;
}

void System::getStatistics()
{
    /*
     * Class instance for sampling statistics from our system.
     */
    double averagedGammaSquared = 0;
    // Performing an average over the Monte Carlo obtained values
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        m_averagedGamma += m_Gamma[alpha];
        averagedGammaSquared += m_Gamma[alpha]*m_Gamma[alpha];
    }
    m_averagedGamma /= double(m_NCf);
    averagedGammaSquared /= double(m_NCf);
    m_varianceGamma = (averagedGammaSquared - m_averagedGamma*m_averagedGamma)/double(m_NCf);
    m_stdGamma = sqrt(m_varianceGamma);
    cout << m_averagedGamma << ", std = " << m_stdGamma << ", variance = " << m_varianceGamma << endl;
}

void System::writeDataToFile(std::string filename)
{
    /*
     * For writing the raw Gamma data to file.
     * Arguments:
     * - filename
     */
    std::ofstream file;
    file.open(m_outputFolder + filename + ".dat");
    file << "acceptanceCounter " << getAcceptanceRate() << endl;
    file << "NCor " << m_NCor << endl;
    file << "NCf " << m_NCf << endl;
    file << "NTherm " << m_NTherm << endl;
    file << "AverageGamma " << m_averagedGamma << endl;
    file << "VarianceGamma " << m_varianceGamma << endl;
    file << "stdGamma " << m_stdGamma << endl;
    if (m_storeThermalizationObservables) {
        for (int i = 0; i < m_NTherm+1; i++) {
            file << m_GammaPreThermalization[i] << endl;
        }
    }
    for (int i = 0; i < m_NCf; i++) {
        file << m_Gamma[i] << endl;
    }
    file.close();
    cout << m_outputFolder + filename + ".dat" << " written" << endl;
}


void System::printAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the metropolis algorithm.
     */
    printf("Acceptancerate: %.16f \n", getAcceptanceRate()); // Times 4 from the Lorentz indices
}

double System::getAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the metropolis algorithm.
     */
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_nUpdates*m_latticeSize*4);
}

void System::writeConfigurationToFile(std::string filename)
{
    /*
     * C-method for writing out configuration to file.
     * Arguments:
     * - filename
     */
    FILE *file; // C method
    file = fopen((m_outputFolder + filename + ".bin").c_str(), "wb");
    for (int t = 0; t < m_N_T; t++) {
        for (int z = 0; z < m_N; z++) {
            for (int y = 0; y < m_N; y++) {
                for (int x = 0; x < m_N; x++) {
                    fwrite(&m_lattice[index(x,y,z,t,m_N,m_N_T)],sizeof(Links),1,file);
                }
            }
        }
    }
    fclose(file);
    cout << m_outputFolder + filename + ".bin" + " written" << endl;
}

void System::loadFieldConfiguration(std::string filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */
    FILE *file; // C method
    file = fopen((m_outputFolder + filename + ".bin").c_str(), "rb"); // CAREFULL HERE!

    for (int t = 0; t < m_N_T; t++) {
        for (int z = 0; z < m_N; z++) {
            for (int y = 0; y < m_N; y++) {
                for (int x = 0; x < m_N; x++) {
//                    cout << "Index: " << x << " " << y << " " << z << " " << t << " Mem.pos.: " << index(x,y,z,t,m_N,m_N_T) << endl;
                    fread(&m_lattice[index(x,y,z,t,m_N,m_N_T)],sizeof(Links),1,file);
                    if (feof(file)) {
                        cout << "EOF!" << endl;
                        exit(0);
                    }
                }
            }
        }
    }

    fclose(file);
    cout << m_outputFolder + filename + ".bin" + " loaded" << endl;
}
