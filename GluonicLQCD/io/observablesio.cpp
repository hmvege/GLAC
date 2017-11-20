#include "observablesio.h"

#include <fstream>
#include <iostream>
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
                          + "observables/"
                          + Parameters::getBatchName() + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "acceptanceCounter " << acceptanceRate << endl;
        file << "NCor " << Parameters::getNCor() << endl;
        file << "NCf " << Parameters::getNCf() << endl;
        file << "NTherm " << Parameters::getNTherm() << endl;
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
        double a = Parameters::getLatticeSpacing();
        double flowStep = Parameters::getFlowEpsilon();
        std::ofstream file;
        std::string fname = Parameters::getFilePath()
                          + Parameters::getOutputFolder()
                          + "flow_observables/"
                          + Parameters::getBatchName() + "_"
                          + observableName + "_flow_config" + std::to_string(configNumber) + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "NFlows " << Parameters::getNFlows() << endl;
        file << "FlowEpsilon " << Parameters::getFlowEpsilon() << endl;
        for (int i = 0; i < Parameters::getNFlows(); i++) {
            file << i << " " << a*sqrt(8*flowStep*(i+1)) << " " << observables[i] << endl;
        }
        file.close();
        if (Parameters::getVerbose()) printf("\n%s written.",fname.c_str());
        std::setprecision(oldPrecision);
    }
}
