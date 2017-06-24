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

using std::cout;
using std::endl;

Metropolis::Metropolis(int N, int N_T, int NCf, int NCor, int Therm, double a, double L, double seed, Correlator *correlator, Action *S)
{
    /*
     * Class for calculating correlators using the Metropolis algorithm.
     * Takes an action object as well as a Gamma functional to be used in the action.
     */
    m_N = N; // Spatial dimensions
    m_N_T = N_T; // Time dimensions
    m_latticeSize = m_N*m_N*m_N*m_N_T;
    m_NCf = NCf; // Number of configurations to run for
    m_NCor = NCor;
    m_NTherm = Therm;
    m_a = a;
    m_L = L;
    setAction(S);
    setCorrelator(correlator);
    m_lattice = new Links[m_latticeSize]; // Lattice, contigious memory allocation
    Gamma = new double[m_NCf]; // Correlator values
    GammaSquared = new double[m_NCf];

    std::mt19937_64 gen(seed); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uni_dist(0,1);
    m_generator = gen;
    m_uniform_distribution = uni_dist;
//    GammaVariance = new double[NCf];
//    GammaStd= new double[NCf];
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        Gamma[alpha] = 0;
        GammaSquared[alpha] = 0;
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
    delete [] Gamma;
    delete [] GammaSquared;
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
//            lattice[i].U[mu] = m_SU3Generator->generate();
            m_lattice[i].U[mu] = m_SU3Generator->updateMatrix();
//            lattice[i].U[mu] = m_SU3Generator->generateIdentity(); // GENERATES IDENTITY FOR TEST! ONE OBSERVABLE SHOULD EQUAL 1!!
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
//    SU3 X = m_SU3Generator->generate(); // Generates a random matrix, SHOULD BE MODIFIED TO X = RST, page 83 Gattinger & Lang
    SU3 X = m_SU3Generator->updateMatrix();
    m_updatedMatrix = X*m_lattice[latticeIndex].U[mu];
}

void Metropolis::update()
{
    /*
     * Sweeps the entire Lattice, and gives every matrix a chance to update.
     */
    for (int i = 0; i < m_N; i++) {
        for (int j = 0; j < m_N; j++) {
            for (int k = 0; k < m_N; k++) {
                for (int l = 0; l < m_N_T; l++) {
                    for (int mu = 0; mu < 4; mu++) {
                        m_S->computeStaple(m_lattice,i,j,k,l,mu);
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(index(i, j, k, l, m_N), mu);
                            m_deltaS = m_S->getDeltaAction(m_lattice, m_updatedMatrix, i, j, k, l, mu);
                            m_expDeltaS = exp(-m_deltaS);
                            if (m_uniform_distribution(m_generator) <= m_expDeltaS)
                            {
                                m_lattice[index(i, j, k, l, m_N)].U[mu].copy(m_updatedMatrix);
                            }
                            else
                            {
                                acceptanceCounter++;
                            }
                        }
                    }
                }
            }
        }
    }
}


void Metropolis::runMetropolis()
{
    cout << "Pre-thermialization correlator: " << m_correlator->calculate(m_lattice) << endl;
    writeConfigurationToFile();
//    exit(1);
    // Running thermalization
    for (int i = 0; i < m_NTherm * m_NCor; i++)
    {
        update();
    }
    cout << "Post-thermialization correlator: " << m_correlator->calculate(m_lattice) << endl;
    cout << "Termalization complete. Acceptance rate: " << getAcceptanceRate() << endl;
//    exit(1);
    // Setting the Metropolis acceptance counter to 0 in order not to count the thermalization
    acceptanceCounter = 0;
    // Main part of algorithm
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        for (int i = 0; i < m_NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
        }
        Gamma[alpha] = m_correlator->calculate(m_lattice);
//        cout << Gamma[alpha] <<endl;|
    }
    cout << "Metropolis completed, line 126" << endl;
    cout << m_correlator->calculate(m_lattice) << endl;
    writeConfigurationToFile();
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

    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        GammaSquared[alpha] = Gamma[alpha]*Gamma[alpha];
    }

    // Performing an average over the Monte Carlo obtained values
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        averagedGamma += Gamma[alpha];
        varianceGamma += (GammaSquared[alpha] - Gamma[alpha]*Gamma[alpha])/double(m_NCf);
    }
    averagedGamma  /= double(m_NCf);
    varianceGamma  /= double(m_NCf);
    stdGamma = sqrt(varianceGamma);
    cout << averagedGamma << " +/- " << stdGamma << " " << varianceGamma << endl;
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
    file << "AverageGamma " << averagedGamma << endl;
    file << "VarianceGamma " << varianceGamma << endl;
    file << "stdGamma " << stdGamma << endl;
    for (int alpha = 0; alpha < m_NCf; alpha++)
    {
        file << Gamma[alpha] << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}


void Metropolis::printAcceptanceRate()
{
    printf("Acceptancerate: %.16f \n", getAcceptanceRate()); // Times 4 from the Lorentz indices
}

void Metropolis::writeConfigurationToFile()
{
    std::ofstream file;
    file.open("../output/configs.txt");
    for (int i = 0; i < m_N; i++) {
        for (int j = 0; j < m_N; j++) {
            for (int k = 0; k < m_N; k++) {
                for (int l = 0; l < m_N_T; l++) {
                    for (int mu = 0; mu < 4; mu++) {
                        for (int i_x = 0; i_x < 9; i_x++) {
                            file << m_lattice[index(i,j,k,l,m_N)].U[mu].mat[i_x].re << " " << m_lattice[index(i,j,k,l,m_N)].U[mu].mat[i_x].im << "      ";
                        }
                        file << endl;
                    }
                    file << endl;
                }
                file << endl;
            }
            file << endl;
        }
        file << endl;
    }
    file.close();
    cout << "configs.txt written" << endl;
}

double Metropolis::getAcceptanceRate()
{
    return double(acceptanceCounter)/double(m_NCf*m_NCor*m_nUpdates*m_latticeSize*4);
}

void Metropolis::loadFieldConfiguration(const char *filename)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     */
    std::ifstream file;
    file.open(filename);
    file.close();
}
