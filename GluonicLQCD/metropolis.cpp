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

Metropolis::Metropolis(int new_N, int new_NCf, int new_NCor, int new_Therm, double new_a)
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
}
Metropolis::~Metropolis()
{
    /*
     * Class destructor
     */
    delete [] lattice;
    delete [] linkMatrices;
    for (int i = 0; i < NCf; i++) { delete [] Gamma[i]; }
    delete [] Gamma;
    delete [] averagedGamma;
    delete [] averagedGammaSquared;
    delete [] varianceGamma;
    delete [] stdGamma;
    delete [] deltaE;
    delete [] deltaE_std;
}

//arma::mat Metropolis::generateMatrix(std::mt19937_64 &gen, std::uniform_real_distribution<double> &randDistr)
//{
//    /*
//     * Most likely try to place inside a class later
//     */
//    int colSize = 3;
//    // Generate a column vector with random numbers
//    double * col1 = new double[colSize];
//    double * col2 = new double[colSize];
//    double * col3 = new double[colSize];
//    for (int i = 0; i < colSize; i++)
//    {
//        col1[i] = randDistr(gen);
//        col2[i] = randDistr(gen);
//    }
//    // Normalize
//    normalizeVector(col1,colSize);
//    // Gram-Schmitt for first column
//    gramSchmitt(col1,col2,colSize);
//    double tempSum = 0;
//    for (int i = 0; i < colSize; i++) { tempSum += col1[i]*col2[i]; }
//    cout << "tempSum " << tempSum << endl;
//    exit(1);
//    delete [] col1;
//    delete [] col2;
//    delete [] col3;
//}

//void Metropolis::normalizeVector(double * v, int n)
//{
//    /*
//     * Normalizes vectors of length n
//     */
//    double vSum = 0;
//    for (int i = 0; i < n; i++) { vSum += v[i]; }
//    for (int i = 0; i < n; i++) { v[i] /= vSum; }
//}


//void Metropolis::testOrthogonality(arma::mat M)
//{
//    cout << M.t()*M << endl;
//}

void Metropolis::latticeSetup(SU3MatrixGenerator *SU3Generator)
{
    /*
     * Sets up the lattice and its matrices.
     */
    lattice = new double[latticeSize];
    linkMatrices = new Links[latticeSize];

//    std::mt19937_64 generator(std::time(nullptr)); // Starting up the Mersenne-Twister19937 function
//    std::uniform_real_distribution<double> distribution(-1,1);

    for (int i = 0; i < latticeSize; i++)
    {
        lattice[i] = 0;
        linkMatrices[i].M0 = SU3Generator->generate();
        linkMatrices[i].M1 = SU3Generator->generate();
        linkMatrices[i].M2 = SU3Generator->generate();
        linkMatrices[i].M3 = SU3Generator->generate();
    }
}

void Metropolis::update(double *x,
                        std::mt19937_64 &gen,
                        std::uniform_real_distribution<double> &epsilon_distribution,
                        std::uniform_real_distribution<double> &uniform_distribution)
{
    /*
     * Private function used for updating our system. Performs the Metropolis algorithm
     */
    for (int i = 0; i < N; i++)
    {
        double x_prev = x[i];
        double oldS = S->getAction(x, i);
        x[i] += epsilon_distribution(gen); // setting a new possible x-position to test for
        double deltaS = S->getAction(x, i) - oldS;
        if ((deltaS > 0) && (exp(-deltaS) < uniform_distribution(gen)))
        {
            x[i] = x_prev;
        }
        else
        {
            acceptanceCounter++;
        }
    }
}

void Metropolis::runMetropolis()
{
    // Setting up random generators
    std::mt19937_64 generator(std::time(nullptr)); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> epsilon_distribution(-epsilon, epsilon);
    std::uniform_real_distribution<double> uniform_distribution(0,1);

    // Setting up array for Gamma-functional values
    Gamma = new double*[NCf];
    for (int i = 0; i < NCf; i++) { Gamma[i] = new double[N]; }
    for (int i = 0; i < NCf; i++) { for (int j = 0; j < N; j++) { Gamma[i][j] = 0; } } // Setting matrix elements to zero

    // Setting up array
    double * x = new double[N]; // Only need one array, as it will always be updated. Note, it is 1D
    for (int i = 0; i < N; i++) { x[i] = 0; }

    // Running thermalization
    for (int i = 0; i < NTherm * NCor; i++)
    {
        update(x, generator, epsilon_distribution, uniform_distribution);
    }

    // Setting the Metropolis acceptance counter to 0 in order not to count the thermalization
    acceptanceCounter = 0;

    // Main part of algorithm
    for (int alpha = 0; alpha < NCf; alpha++)
    {
        for (int i = 0; i < NCor; i++) // Updating NCor times before updating the Gamma function
        {
            update(x, generator, epsilon_distribution, uniform_distribution);
        }
        for (int n = 0; n < N; n++)
        {
            Gamma[alpha][n] = gammaFunctional(x,n,N);
        }
    }

    // De-allocating arrays
    delete [] x;
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
