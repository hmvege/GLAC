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
                                               double *observables,
                                               double *thermalizationObservables,
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

void IO::writeFlowObservableToFile(double *observableStatistics,
                                                  double *t,
                                                  double *observables,
                                                  std::string observableName)
{
    /*
     * Method for writing a flow variable to file.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(15);
        std::ofstream file;
        std::string fname = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName() + "_" + observableName + "_flow.dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "NFlows " << Parameters::getNFlows() << endl;
        file << "FlowEpsilon " << Parameters::getFlowEpsilon() << endl;
        file << "Average" << observableName << " " << observableStatistics[0] << endl;
        file << "Variance" << observableName << " " << observableStatistics[1] << endl;
        file << "std" << observableName << " " << observableStatistics[2] << endl;
        for (int i = 0; i < Parameters::getNFlows(); i++) {
            file << i << " " << t[i] << " " << observables[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
        std::setprecision(oldPrecision);
    }
}
