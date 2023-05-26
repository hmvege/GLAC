#include "observablesio.h"
#include "config/parameters.h"
#include "parallelization/communicator.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

/*!
 * \brief IO::writeObservablesToFile method for writing an observable to file, using .15 precision.
 * \param acceptanceRate acceptance rate from the Metropolis algorithm.
 * \param averagedObservable the average of the observable.
 * \param varianceObservable the variance of the observable.
 * \param stdObservable the standard deviation of the observable.
 * \param observables vector of doubles containing the observable.
 * \param NObs the number of configurations, i.e. observables, to be written to file.
 * \param observableName name of the observable. Topological charge, plaquette, energy ect.
 *
 * The top 8 rows contains,
 * - beta,
 * - acceptance rate,
 * - number of correlation updates used,
 * - number of configurations,
 * - number of thermalization steps,
 * - the average of the observable,
 * - the variance of the observable,
 * - the standard deviation of the observable,
 */
void IO::writeObservablesToFile(const double acceptanceRate,
                                const double averagedObservable,
                                const double varianceObservable,
                                const double stdObservable,
                                const std::vector<double> &observables,
                                const unsigned int NObs,
                                const std::string &observableName)
{
    /*
     * Function to write thermalization and configuration observables to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0)
    {
        auto oldPrecision = cout.precision(15);
        std::ofstream file;
        std::string fname = Parameters::getFilePath()
                          + Parameters::getOutputFolder()
                          + Parameters::getBatchName() + "/"
                          + "observables/" + observableName + "_"
                          + Parameters::getBatchName() + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "acceptanceCounter " << acceptanceRate << endl;
        file << "NCor " << Parameters::getNCor() << endl;
        file << "NCf " << Parameters::getNCf() << endl;
        file << "NTherm " << Parameters::getNTherm() << endl;
        file << std::fixed << std::setprecision(15);
        file << "Average" << observableName << " " << averagedObservable << endl;
        file << "Variance" << observableName << " " << varianceObservable << endl;
        file << "std"  << observableName << " " << stdObservable << endl;
        for (unsigned int i = 0; i < NObs; i++){
            file << observables[i] << endl;
        }
        file.close();
        printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}

/*!
 * \brief IO::writeFlowObservableToFile method for writing a flow observable to file, using .15 precision.
 * \param observables a vector of observables.
 * \param observableName observable name.
 * \param configNumber the configuration number.
 *
 * Top 3 rows of the output file will contain,
 * - beta,
 * - NFlows,
 * - FlowEpsilon,
 * with the rest containing the flowed observable.
 */
void IO::writeFlowObservableToFile(const std::vector<double> &observables,
                                   const std::string &observableName,
                                   const unsigned int configNumber)
{
    /*
     * Method for writing a flow variable to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0)
    {
        auto oldPrecision = cout.precision(15);
        double flowStep = Parameters::getFlowEpsilon();
        std::ofstream file;

        // Converting config number to a more machine friendly layout
        char cfg_number[6];
        sprintf(cfg_number,"%05d",configNumber + Parameters::getConfigStartNumber());

        // Sets up the file name
        std::string fname = Parameters::getOutputFolder()
                          + Parameters::getBatchName() + "/"
                          + "flow_observables/" + observableName + "/"
                          + Parameters::getBatchName() + "_"
                          + observableName + "_flow_config" + std::string(cfg_number) + ".dat";

        // Opens file
        file.open(Parameters::getFilePath() + fname);

        // Writes out essential information about the observable
        file << "beta " << Parameters::getBeta() << endl;
        file << "NFlows " << Parameters::getNFlows() << endl;
        file << "FlowEpsilon " << Parameters::getFlowEpsilon() << endl;

        // Sets temporary high precision
        file << std::fixed << std::setprecision(15);

        // Writes out the observable
        for (unsigned int i = 0; i < Parameters::getNFlows() + 1; i++)
        {
            file << i*flowStep << " " << observables[i] << endl;
        }

        // Closes file, prints file written if verbose and reverts the precision
        file.close();
        if (Parameters::getVerbose()) printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}

/*!
 * \brief IO::writeMatrixToFile writes out a matrix of observables to file, i.e. the topological charge in Euclidean time. Uses a .15 precision.
 * \param observables the observables in a vector of doubles.
 * \param observableName the name of the observable.
 * \param configNumber the configuration number.
 * \param N the number of columns in the matrix of observables. The number of rows in the matrix is the number of flow steps, NFlows, retrieved from Parameters::getNFlows().
 *
 * Top 3 rows of the output file will contain,
 * - beta,
 * - NFlows,
 * - FlowEpsilon,
 * with the rest containing the flowed observable.
 */
void IO::writeMatrixToFile(const std::vector<double> &observables, const std::string &observableName,
                           const unsigned int configNumber, const unsigned int N)
{
    /*
     * Function for writing out a matrix of observables to file, e.g. the topc in euclidean time.
     */

    // Writes to file
    if (Parallel::Communicator::getProcessRank() == 0)
    {
        auto oldPrecision = cout.precision(15);
        double flowStep = Parameters::getFlowEpsilon();
        std::ofstream file;

        // Converting config number to a more machine friendly layout
        char cfg_number[6];
        sprintf(cfg_number,"%05d",configNumber + Parameters::getConfigStartNumber());

        // Sets up the file name
        std::string fname = Parameters::getOutputFolder()
                          + Parameters::getBatchName() + "/"
                          + "flow_observables/" + observableName + "/"
                          + Parameters::getBatchName() + "_"
                          + observableName + "_flow_config" + std::string(cfg_number) + ".dat";

        // Opens file
        file.open(Parameters::getFilePath() + fname);

        // Writes out essential information about the observable
        file << "beta " << Parameters::getBeta() << endl;
        file << "NFlows " << Parameters::getNFlows() << endl;
        file << "FlowEpsilon " << Parameters::getFlowEpsilon() << endl;

        // Sets temporary high precision
        file << std::fixed << std::setprecision(15);

        // Writing out the matrix columns along each line
        for (unsigned int iFlow = 0; iFlow < Parameters::getNFlows() + 1; iFlow++)
        {
            file << iFlow*flowStep;
            for (unsigned int i = 0; i < N; i++)
            {
                file << " " << observables[iFlow * N + i];
            }
            file << endl;
        }

        // Closes file, prints file written if verbose and reverts the precision
        file.close();
        if (Parameters::getVerbose()) printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}
