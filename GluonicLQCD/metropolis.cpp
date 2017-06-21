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
#include "su3matrixgenerator.h"

using std::cout;
using std::endl;

Metropolis::Metropolis(int new_N, int new_NCf, int new_NCor, int new_Therm, double new_a, double new_L, Correlator *new_correlator, Action *new_S)
{
    /*
     * Class for calculating correlators using the Metropolis algorithm.
     * Takes an action object as well as a Gamma functional to be used in the action.
     */
    N = new_N; // Lattice sites
    latticeSize = N*N*N*N;
    NCf = new_NCf; // Number of configurations to run for
    NCor = new_NCor;
    NTherm = new_Therm;
    a = new_a;
    L = new_L;
    setAction(new_S);
    setCorrelator(new_correlator);
    lattice = new Links[latticeSize]; // Lattice, contigious memory allocation
    Gamma = new double[NCf]; // Correlator values
    GammaSquared = new double[NCf];
//    GammaVariance = new double[NCf];
//    GammaStd= new double[NCf];
    for (int alpha = 0; alpha < NCf; alpha++)
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
    delete [] lattice;
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
    for (int i = 0; i < latticeSize; i++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            lattice[i].U[mu] = m_SU3Generator->generate();
//            lattice[i].U[mu] = m_SU3Generator->generateIdentity(); // GENERATES IDENTITY FOR TEST! ONE OBSERVABLE SHOULD EQUAL 1!!
        }
    }
}

void Metropolis::updateLink(int i, int mu)
{
    /*
     * Private function used for updating our system. Updates a single gauge link.
     * Arguments:
     *  i   : spacetime index
     *  mu  : Lorentz index
     */
    SU3 X = m_SU3Generator->generate(); // Generates a random matrix, SHOULD BE MODIFIED TO X = RST, page 83 Gattinger & Lang
    updatedMatrix = X*lattice[i].U[mu];
}

void Metropolis::update()
{
    // Updates the entire Lattice
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < N; l++) {
                    for (int mu = 0; mu < 4; mu++)
                    {
                        S->computeStaple(lattice,i,j,k,l,mu);
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(index(i, j, k, l, N), mu);
                            deltaS = S->getDeltaAction(lattice, updatedMatrix, i, j, k, l, mu);
                            expDeltaS = exp(-deltaS);
                            if (m_uniform_distribution(m_generator) <= expDeltaS)
                            {
                                lattice[i].U[mu].copy(updatedMatrix);
                            }
                            else
                            {
                                acceptanceCounter++;
                            }
                        }
//                        cout << m_correlator->calculate(lattice) << endl;
                        // Calculate action
//                        cout << deltaS << endl;
//                        if ((expDeltaS >= 1) || (m_uniform_distribution(m_generator) <= expDeltaS))
                    }
                }
            }
        }
    }
}


void Metropolis::runMetropolis()
{
    cout << "Running Metropolis for Gluon action. Line 107" << endl;
    // Initializing storage variables

    // Running thermalization

//    cout << m_correlator->calculate(lattice) << endl;
    for (int i = 0; i < NTherm * NCor; i++)
    {
        update();
    }
    cout << "Termalization complete. Line 116." << endl;
//    printAcceptanceRate();
//    exit(1);
//    exit(1);
    // Setting the Metropolis acceptance counter to 0 in order not to count the thermalization
    acceptanceCounter = 0;

    // Main part of algorithm
    for (int alpha = 0; alpha < NCf; alpha++)
    {
        for (int i = 0; i < NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
        }
        Gamma[alpha] = m_correlator->calculate(lattice);
    }
    cout << "Metropolis completed, line 126" << endl;
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

    for (int alpha = 0; alpha < NCf; alpha++)
    {
        GammaSquared[alpha] = Gamma[alpha]*Gamma[alpha];
    }

    // Performing an average over the Monte Carlo obtained values
    for (int alpha = 0; alpha < NCf; alpha++)
    {
        averagedGamma += Gamma[alpha];
        varianceGamma += (GammaSquared[alpha] - Gamma[alpha]*Gamma[alpha])/double(NCf);
    }
    averagedGamma  /= double(NCf);
    varianceGamma  /= double(NCf);
    stdGamma = sqrt(varianceGamma);
    cout << averagedGamma << " +/- " << stdGamma << endl;
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
    file << "acceptanceCounter " << double(acceptanceCounter)/double(NCf*NCor*latticeSize*4) << endl;
    file << "NCor " << NCor << endl;
    file << "NCf " << NCf << endl;
    file << "NTherm " << NTherm << endl;
    file << "AverageGamma " << averagedGamma << endl;
    file << "VarianceGamma " << varianceGamma << endl;
    file << "stdGamma " << stdGamma << endl;
    for (int alpha = 0; alpha < NCf; alpha++)
    {
        file << Gamma[alpha] << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}

//void Metropolis::writeStatisticsToFile(const char *filename)//, double * dE, double * averagedGamma, double * averagedGammaSquared, int acceptanceCounter)
//{
    /*
     * Writes statistics to file about:
     * acceptanceCounter:   the number of accepted configurations
     * NCor:                number of times between each sampling of the functional
     * NCf:                 number of paths we are looking at
     * t=n*a:               points on lattice
     * dE:                  energy for a given point on lattice
     * dE error:            the error in the dE measurment
     * variance:            Var(G)
     * standardDeviation:   std(G)
     */
//    std::ofstream file;
//    file.open(filename);
//    file << "acceptanceCounter " << double(acceptanceCounter)/double(NCf*NCor*latticeSize*4) << endl;
//    file << "NCor " << NCor << endl;
//    file << "NCf " << NCf << endl;
//    file << "NTherm " << NTherm << endl;
//    file.close();
//    cout << filename << " written" << endl;
//}

//void Metropolis::printEnergies()
//{
//    /*
//     * Printing the energies in the calculation
//     */
//    for (int n = 0; n < N; n++)
//    {
//        cout << deltaE[n] << endl;
//    }
//}

void Metropolis::printAcceptanceRate()
{
    printf("Acceptancerate: %.16f \n", double(acceptanceCounter)/double(NCf*NCor*latticeSize*4)); // Times 4 from the Lorentz indices
}


void Metropolis::writeConfigurationToFile()
{
    std::ofstream file;
    file.open("../output/configs.txt");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < N; l++) {
                    for (int mu = 0; mu < 4; mu++) {
                        for (int i_x = 0; i_x < 9; i_x++) {
                            file << lattice[index(i,j,k,l,N)].U[mu].mat[i_x] << " ";
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
