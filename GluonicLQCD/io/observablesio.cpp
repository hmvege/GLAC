#include "observablesio.h"
#include "config/parameters.h"
#include "parallelization/communicator.h"
#include <cmath>
#include <fstream>
#include <iomanip>

void IO::writeObservablesToFile(double acceptanceRate,
                                double averagedObservable,
                                double varianceObservable,
                                double stdObservable,
                                double *observables,
                                int NObs,
                                std::string observableName)
{
    /*
     * Function to write thermalization and configuration observables to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
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
        for (int i = 0; i < NObs; i++) {
            file << observables[i] << endl;
        }
        file.close();
        printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}

void IO::writeFlowObservableToFile(double *observables,
                                   std::string observableName,
                                   int configNumber)
{
    /*
     * Method for writing a flow variable to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
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
        for (int i = 0; i < Parameters::getNFlows() + 1; i++) {
            file << i*flowStep << " " << observables[i] << endl;
        }

        // Closes file, prints file written if verbose and reverts the precision
        file.close();
        if (Parameters::getVerbose()) printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}

void IO::writeMatrixToFile(double * observables, std::string observableName, int configNumber, int N)
{
    /*
     * Function for writing out a matrix of observables to file, e.g. the topc in euclidean time.
     */

    // Writes to file
    if (Parallel::Communicator::getProcessRank() == 0) {
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


//        printf("\nMember to remove memory check in observables io line 127\n");

        // Writing out the matrix columns along each line
        for (int iFlow = 0; iFlow < Parameters::getNFlows() + 1; iFlow++) {
            file << iFlow*flowStep;
            for (int i = 0; i < N; i++) {
                file << " " << observables[iFlow * N + i];

//                // MEMORY CHECKING
//                if ((iFlow * N + i) > ((Parameters::getNFlows()+1)*Parameters::getNTemporal())) {
//                    cout<< "ERROR! EXCEEDING BOUNDS! N="<< N <<", i="<<i<<", iFlow="<<iFlow<< ": " << iFlow * N + i << " " << Parameters::getNFlows()*Parameters::getNTemporal() + 1 << endl;
//                    exit(1);
//                }
            }
            file << endl;
        }

//        // TEMPTEMP
//        m_tempTopctArray = new double[Parameters::getNTemporal() * Parameters::getNFlows()];
//        for (int iFlow = 0; iFlow < Parameters::getNFlows(); iFlow++) {
//            for (int it = 0; it < Parameters::getNTemporal(); it++) {
//                m_tempTopctArray[iFlow*Parameters::getNTemporal() + it] = 0;
//            }
//        }


//        printf("\nOK AFTER IO WRITING!!\n");

        // Closes file, prints file written if verbose and reverts the precision
        file.close();
        if (Parameters::getVerbose()) printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}
