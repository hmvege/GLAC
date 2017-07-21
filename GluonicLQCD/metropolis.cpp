#include <random>   // For Mersenne-Twister19937
#include <ctime>    // For random seed
#include <cmath>    // For exp()
#include <fstream>
#include <iostream>
#include "metropolis.h"
#include "actions/action.h"
#include "correlators/correlator.h"
#include "functions.h"
#include "links.h"
#include "matrices/su3matrixgenerator.h"

//TEMP
#include "unittests.h"

#include <cstdio>
#include <cstdlib>

using std::cout;
using std::endl;

Metropolis::Metropolis(int N, int N_T, int NCf, int NCor, int NTherm, double a, double L, double seed, Correlator *correlator, Action *S)
{
    /*
     * Class for calculating correlators using the Metropolis algorithm.
     * Takes an action object as well as a Gamma functional to be used in the action.
     */
    m_N = N; // Spatial dimensions
    m_N_T = N_T; // Time dimensions
    m_latticeSize = N*N*N*N_T;
    m_NCf = NCf; // Number of configurations to run for
    m_NCor = NCor;
    m_NTherm = NTherm;
    m_a = a;
    m_L = L;
    setAction(S);
    setCorrelator(correlator);
    m_lattice = new Links[m_latticeSize]; // Lattice, contigious memory allocation
    m_Gamma = new double[m_NCf]; // Correlator values
    m_GammaSquared = new double[m_NCf];

    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uni_dist(0,1);
    m_generator = gen;
    m_uniform_distribution = uni_dist;
//    GammaVariance = new double[NCf];
//    GammaStd= new double[NCf];
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        m_Gamma[alpha] = 0;
        m_GammaSquared[alpha] = 0;
//        GammaVariance[alpha] = 0;
//        GammaStd[alpha] = 0;
    }
}

Metropolis::~Metropolis()
{
    /*
     * Class destructor
     */
    delete [] m_lattice;
    delete [] m_Gamma;
    delete [] m_GammaSquared;
//    delete [] GammaVariance;
//    delete [] GammaStd;
}

void Metropolis::latticeSetup(SU3MatrixGenerator *SU3Generator)
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
            m_lattice[i].U[mu] = m_SU3Generator->generateIdentity(); // GENERATES IDENTITY FOR TEST! ONE OBSERVABLE SHOULD EQUAL 1!!
        }
    }
}

void Metropolis::updateLink(int latticeIndex, int mu)
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

void Metropolis::update()
{
    /*
     * Sweeps the entire Lattice, and gives every matrix a chance to update.
     */
//    double updateS = 0; //TEMP
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
//                            cout << "Do I forget to divide by 2 somewhere?" << endl;
//                            updateS+=m_deltaS; // For checking if something is wrong in the update
                            if (exp(-m_deltaS) > m_uniform_distribution(m_generator))
                            {
                                m_lattice[index(x, y, z, t, m_N, m_N_T)].U[mu].copy(m_updatedMatrix);
                            }
                            else
                            {
                                m_acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
//    cout << updateS/double(4*m_latticeSize*m_nUpdates) << endl;
}


void Metropolis::runMetropolis()
{
//    loadFieldConfiguration("conf0.bin");
    cout << "Pre-thermialization correlator:  " << m_correlator->calculate(m_lattice) << endl;
    // Running thermalization
    for (int i = 0; i < m_NTherm * m_NCor; i++)
    {
        update();
        // Print correlator every somehting
        if ((i+1) % 10 == 0) {
            cout << m_correlator->calculate(m_lattice) << endl;
        }
    }
    cout << "Post-thermialization correlator: " << m_correlator->calculate(m_lattice) << endl;
    cout << "Termalization complete. Acceptance rate: " << m_acceptanceCounter/double(4*m_latticeSize*m_nUpdates*m_NTherm*m_NCor) << endl;
    // Setting the Metropolis acceptance counter to 0 in order not to count the thermalization
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
    cout << m_correlator->calculate(m_lattice) << endl;
}

void Metropolis::sampleSystem()
{
    /*
     * For sampling statistics/getting correlators
     */
    cout << "Not implemented yet." << endl;
}

void Metropolis::getStatistics()
{
    /*
     * Class instance for sampling statistics from our system.
     */
//    deltaE          = new double[N];
//    deltaE_std      = new double[N];
    double averagedGammaSquared = 0;
//    for (int alpha = 0; alpha < m_NCf; alpha++)
//    {
//        GammaSquared[alpha] = Gamma[alpha]*Gamma[alpha];
//    }
    // Performing an average over the Monte Carlo obtained values
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        m_averagedGamma += m_Gamma[alpha];
        averagedGammaSquared += m_Gamma[alpha]*m_Gamma[alpha];
    }
//    for (int alpha = 0; alpha < m_NCf; alpha++)
//    {
//        varianceGamma += (GammaSquared[alpha] - Gamma[alpha]*Gamma[alpha])/double(m_NCf);
//    }
    m_averagedGamma /= double(m_NCf);
    averagedGammaSquared /= double(m_NCf);
    m_varianceGamma = (averagedGammaSquared - m_averagedGamma*m_averagedGamma)/double(m_NCf);
    m_stdGamma = sqrt(m_varianceGamma);
    cout << m_averagedGamma << ", std = " << m_stdGamma << ", variance = " << m_varianceGamma << endl;
    // Getting change in energy & calculating variance & standard deviation of G
//    for (int n = 0; n < N; n++)
//    {
//        deltaE[n] = log(averagedGamma[n]/averagedGamma[(n+1) % N])/a;
//    }
//    averagedGamma[n] /= double(NCf);

    // Calculating the uncertainty in dE(hand calculation for analytic expression done beforehand)
//    for (int n = 0; n < N; n++)
//    {
//        deltaE_std[n] = sqrt(pow(stdGamma[n]/averagedGamma[n],2) + pow(stdGamma[(n+1)%N]/averagedGamma[(n+1)%N],2))/a;
//    }
//    delete [] deltaE;
//    delete [] deltaE_std;
}


void Metropolis::writeDataToFile(const char *filename)
{
    /*
     * For writing the raw Gamma data to file.
     */
    std::ofstream file;
    file.open(filename);
    file << "acceptanceCounter " << getAcceptanceRate() << endl;
    file << "NCor " << m_NCor << endl;
    file << "NCf " << m_NCf << endl;
    file << "NTherm " << m_NTherm << endl;
    file << "AverageGamma " << m_averagedGamma << endl;
    file << "VarianceGamma " << m_varianceGamma << endl;
    file << "stdGamma " << m_stdGamma << endl;
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        file << m_Gamma[alpha] << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}


void Metropolis::printAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the Metropolis algorithm.
     */
    printf("Acceptancerate: %.16f \n", getAcceptanceRate()); // Times 4 from the Lorentz indices
}

void Metropolis::writeConfigurationToFile(std::string filename)
{
//    std::ofstream file("../output/" + filename, std::ofstream::binary); // Old c++ method
//    for (int t = 0; t < m_N_T; t++) {
//        for (int z = 0; z < m_N; z++) {
//            for (int y = 0; y < m_N; y++) {
//                for (int x = 0; x < m_N; x++) {
//                    for (int mu = 0; mu < 4; mu++) {
//                        file.write(reinterpret_cast<const char*> (&m_lattice[index(x,y,z,t,m_N,m_N_T)].U[mu]), sizeof(SU3));
//                    }
//                    file << endl;
//                }
//                file << endl;
//            }
//            file << endl;
//        }
//        file << endl;
//    }
//    file.close();
    FILE *file; // C method
    file = fopen((m_outputFolder+filename).c_str(), "wb");
    for (int t = 0; t < m_N_T; t++) {
        for (int z = 0; z < m_N; z++) {
            for (int y = 0; y < m_N; y++) {
                for (int x = 0; x < m_N; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fwrite(&m_lattice[index(x,y,z,t,m_N,m_N_T)].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }
    fclose(file);
    cout << m_outputFolder + filename  + " written" << endl;
}

double Metropolis::getAcceptanceRate()
{
    /*
     * Returns the acceptance ratio of the main run of the Metropolis algorithm.
     */
    return double(m_acceptanceCounter)/double(m_NCf*m_NCor*m_nUpdates*m_latticeSize*4);
}

void Metropolis::loadFieldConfiguration(std::string filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     */
    FILE *file; // C method
    file = fopen((m_inputFolder + filename).c_str(), "rb");
    for (int t = 0; t < m_N_T; t++) {
        for (int z = 0; z < m_N; z++) {
            for (int y = 0; y < m_N; y++) {
                for (int x = 0; x < m_N; x++) {
                    for (int mu = 0; mu < 4; mu++) {
                        fread(&m_lattice[index(x,y,z,t,m_N,m_N_T)].U[mu],sizeof(SU3),1,file);
                    }
                }
            }
        }
    }
    fclose(file);
    cout << m_inputFolder + filename  + " loaded" << endl;
}
