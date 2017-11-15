#include "observablesio.h"

#include <fstream>
#include <iostream>
#include <iomanip>

//IO::ObservablesIO::ObservablesIO()
//{
//    // Add static const std::string class_name for checking what observables we will write to file?
//}

//IO::ObservablesIO::~ObservablesIO()
//{
//}

void IO::writeObservablesToFile(double acceptanceRate,
                                               double averagedObservable,
                                               double varianceObservable,
                                               double stdObservable,
                                               double *observables,
                                               std::string observableName)
{
    /*
     * Function to write configuration observables to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(15);
        std::ofstream file;
        std::string fname = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName() + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "acceptanceCounter " << acceptanceRate << endl;
        file << "NCor " << Parameters::getNCor() << endl;
        file << "NCf " << Parameters::getNCf() << endl;
        file << "NTherm " << Parameters::getNTherm() << endl;
        file << "Average" << observableName << " " << averagedObservable << endl;
        file << "Variance" << observableName << " " << varianceObservable << endl;
        file << "std"  << observableName << " " << stdObservable << endl;
        for (int i = 0; i < Parameters::getNCf(); i++) {
            file << observables[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
        std::setprecision(oldPrecision);
    }
}

void IO::writeObservablesToFile(double acceptanceRate,
                                               double averagedObservable,
                                               double varianceObservable,
                                               double stdObservable,
                                               double *thermalizationObservables,
                                               double *observables,
                                               std::string observableName)
{
    /*
     * Function to write thermalization and configuration observables to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(15);
        std::ofstream file;
        std::string fname = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName() + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "acceptanceCounter " << acceptanceRate << endl;
        file << "NCor " << Parameters::getNCor() << endl;
        file << "NCf " << Parameters::getNCf() << endl;
        file << "NTherm " << Parameters::getNTherm() << endl;
        file << "Average" << observableName << " " << averagedObservable << endl;
        file << "Variance" << observableName << " " << varianceObservable << endl;
        file << "std"  << observableName << " " << stdObservable << endl;
        for (int i = 0; i < Parameters::getNTherm() + 1; i++) {
            file << thermalizationObservables[i] << endl;
        }
        for (int i = 0; i < Parameters::getNCf(); i++) {
            file << observables[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
        std::setprecision(oldPrecision);
    }
}

void IO::writeFlowObservableToFile(double averagedObservable,
                                   double varianceObservable,
                                   double stdObservable,
                                   double *observables,
                                   std::string observableName,
                                   int configNumber)
{
    /*
     * Method for writing a flow variable to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(15);
        double a = Parameters::getLatticeSpacing();
        double flowStep = Parameters::getFlowEpsilon();
        std::ofstream file;
        std::string fname = Parameters::getFilePath()
                          + Parameters::getOutputFolder()
                          + Parameters::getBatchName() + "_"
                          + observableName + "_flow_config" + std::to_string(configNumber) + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "NFlows " << Parameters::getNFlows() << endl;
        file << "FlowEpsilon " << Parameters::getFlowEpsilon() << endl;
        file << "Average" << observableName << " " << averagedObservable << endl;
        file << "Variance" << observableName << " " << varianceObservable << endl;
        file << "std" << observableName << " " << stdObservable << endl;
        for (int i = 0; i < Parameters::getNFlows(); i++) {
            file << i << " " << a*sqrt(8*flowStep*(i+1)) << " " << observables[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
        std::setprecision(oldPrecision);
    }
}

void IO::writeDataToFile(double averagedObservable,
                         double varianceObservable,
                         double stdObservable,
                         double *observables,
                         std::string observableName)
{
    /*
     * Method for writing a general observable to file with statistics.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(15);
        std::ofstream file;
        std::string fname = Parameters::getFilePath() +
                            Parameters::getOutputFolder() +
                            Parameters::getBatchName() + "_configNumber" + "_" + observableName + ".dat";
        file.open(fname);
        file << "Average" << observableName << " " << averagedObservable << endl;
        file << "Variance" << observableName << " " << varianceObservable << endl;
        file << "std" << observableName << " " << stdObservable << endl;
        for (int i = 0; i < Parameters::getNFlows(); i++) {
            file << i << " " << observables[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
        std::setprecision(oldPrecision);
    }
}
