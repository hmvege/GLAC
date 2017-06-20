#include <random>   // For Mersenne-Twister19937
#include <ctime>    // For random seed
#include <cmath>    // For exp()
#include <fstream>
#include <iostream>
#include "metropolis.h"
#include "actions/action.h"
#include "functions.h"
#include "links.h"
#include "su3matrixgenerator.h"

using std::cout;
using std::endl;
void printArray(double *x, int N);

Metropolis::Metropolis(int new_N, int new_NCf, int new_NCor, int new_Therm, double new_a, double new_L)
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
    lattice = new Links[latticeSize];
    // Setting up array for Gamma-functional values
    Gamma = new double*[NCf];
    for (int i = 0; i < NCf; i++) { Gamma[i] = new double[N]; }
    for (int i = 0; i < NCf; i++) { for (int j = 0; j < N; j++) { Gamma[i][j] = 0; } } // Setting matrix elements to zero

}
Metropolis::~Metropolis()
{
    /*
     * Class destructor
     */
    delete [] lattice;
    for (int i = 0; i < NCf; i++) { delete [] Gamma[i]; }
    delete [] Gamma;
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
    SU3 X = m_SU3Generator->generate(); // Generates a random matrix, SHOULD BE MODIVIED TO X = RST, page 83 Gattinger & Lang
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
                        for (int n = 0; n < m_nUpdates; n++) // Runs avg 10 updates on link, as that is less costly than other parts
                        {
                            updateLink(index(i, j, k, l, N), mu);
                        }
                        // Calculate action
                        deltaS = exp(-S->getDeltaAction(lattice, updatedMatrix, i, j, k, l, mu));
                        if ((deltaS >= 1) || (m_uniform_distribution(m_generator) <= deltaS))
                        {
                            lattice[i].U[mu] = updatedMatrix;
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


void Metropolis::runMetropolis()
{
    cout << "Running Metropolis for Gluon action. Line 107" << endl;
    // Initializing storage variables

    // Running thermalization
    for (int i = 0; i < NTherm * NCor; i++)
    {
        update();
    }
    cout << "Termalization complete. Line 116. Acceptance rate: " << double(acceptanceCounter)/double( latticeSize*4*NTherm*NCor ) << endl;
    exit(1);
    // Setting the Metropolis acceptance counter to 0 in order not to count the thermalization
    acceptanceCounter = 0;

    // Main part of algorithm
    for (int alpha = 0; alpha < NCf; alpha++)
    {
        for (int i = 0; i < NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update();
//            double oldS = S->getAction(lattice, i);

//            double deltaS = S->getAction(lattice, i) - oldS;
//            if (uniform_distribution(gen) <= exp(-deltaS))
//            {
//                // Accept new configuration
//            }
//            else
//            {
//                acceptanceCounter++;
//            }
        }
//        for (int n = 0; n < N; n++)
//        {
//            Gamma[alpha][n] = gammaFunctional(lattice,n,N);
//        }
    }
}

void Metropolis::getStatistics()
{
    /*
     * Class instance for sampling statistics from our system.
     */
    averagedGamma           = new double[N];
    averagedGammaSquared    = new double[N];
    varianceGamma           = new double[N];
    stdGamma                = new double[N];
    deltaE                  = new double[N];
    deltaE_std              = new double[N];
    for (int i = 0; i < N; i++)
    {
        averagedGamma[i]        = 0;
        averagedGammaSquared[i] = 0;
        varianceGamma[i]        = 0;
        stdGamma[i]             = 0;
        deltaE[i]               = 0;
        deltaE_std[i]           = 0;
    }

    // Performing an average over the Monte Carlo obtained values
    for (int n = 0; n < N; n++)
    {
        for (int alpha = 0; alpha < NCf; alpha++)
        {
            averagedGamma[n] += Gamma[alpha][n];
            averagedGammaSquared[n] += Gamma[alpha][n]*Gamma[alpha][n];
        }
        averagedGamma[n] /= double(NCf);
        averagedGammaSquared[n] /= double(NCf);
    }

    // Getting change in energy & calculating variance & standard deviation of G
    for (int n = 0; n < N; n++)
    {
        varianceGamma[n] = (averagedGammaSquared[n] - averagedGamma[n]*averagedGamma[n])/NCf;
        stdGamma[n] = sqrt(varianceGamma[n]);
        deltaE[n] = log(averagedGamma[n]/averagedGamma[(n+1) % N])/a;
    }

    // Calculating the uncertainty in dE(hand calculation for analytic expression done beforehand)
    for (int n = 0; n < N; n++)
    {
        deltaE_std[n] = sqrt(pow(stdGamma[n]/averagedGamma[n],2) + pow(stdGamma[(n+1)%N]/averagedGamma[(n+1)%N],2))/a;
    }
    delete [] averagedGamma;
    delete [] averagedGammaSquared;
    delete [] varianceGamma;
    delete [] stdGamma;
    delete [] deltaE;
    delete [] deltaE_std;
}

void Metropolis::writeDataToFile(const char *filename)
{
    /*
     * For writing the raw Gamma data to file.
     */
    std::ofstream file;
    file.open(filename);
    file << "acceptanceCounter " << double(acceptanceCounter)/double(N*NCf*NCor) << endl;
    file << "NCor " << NCor << endl;
    file << "NCf " << NCf << endl;
    file << "NTherm " << NTherm << endl;
    for (int i = 0; i < NCf; i++)
    {
        for (int j = 0; j < N; j++)
        {
            file << Gamma[i][j] << " ";
        }
        file << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}

void Metropolis::writeStatisticsToFile(const char *filename)//, double * dE, double * averagedGamma, double * averagedGammaSquared, int acceptanceCounter)
{
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
    std::ofstream file;
    file.open(filename);
    file << "acceptanceCounter " << double(acceptanceCounter)/double(N*NCf*NCor) << endl;
    file << "NCor " << NCor << endl;
    file << "NCf " << NCf << endl;
    file << "NTherm " << NTherm << endl;
    for (int n = 0; n < N; n++)
    {
        file << n*a << " "
             << deltaE[n] << " "
             << deltaE_std[n] << " "
             << varianceGamma[n] << " "
             << stdGamma[n] << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}

void Metropolis::printEnergies()
{
    /*
     * Printing the energies in the calculation
     */
    for (int n = 0; n < N; n++)
    {
        cout << deltaE[n] << endl;
    }
}

void Metropolis::printAcceptanceRate()
{
    printf("Acceptancerate: %f \n", double(acceptanceCounter)/double(NCf*NCor*N));
}

void printArray(double *x, int N)
{
    for (int i = 0; i < N; i++)
    {
        printf("%.10f \n", x[i]);
    }
}
