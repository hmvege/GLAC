#include "observablesio.h"

#include <fstream>
#include <iostream>
#include <iomanip>

IO::ObservablesIO::ObservablesIO()
{

}

void IO::ObservablesIO::writeObservablesToFile()
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        std::ofstream file;
        std::string fname = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName() + ".dat";
        file.open(fname);
        file << "beta " << Parameters::getBeta() << endl;
        file << "acceptanceCounter " << getAcceptanceRate() << endl;
        file << "NCor " << Parameters::getNCor() << endl;
        file << "NCf " << Parameters::getNCf() << endl;
        file << "NTherm " << Parameters::getNTherm() << endl;
        file << std::setprecision(15) << "AverageObservable " << m_averagedObservable << endl;
        file << std::setprecision(15) << "VarianceObservable " << m_varianceObservable << endl;
        file << std::setprecision(15) << "stdObservable " << m_stdObservable << endl;
        if (m_storeThermalizationObservables) {
            for (int i = 0; i < m_NTherm+1; i++) {
                file << std::setprecision(15) << m_observablePreThermalization[i] << endl;
            }
            file << endl;
        }
        for (int i = 0; i < m_NCf; i++) {
            file << std::setprecision(15) << m_observable[i] << endl;
        }
        file.close();
        cout << fname << " written." << endl;
    }
}

void IO::ObservablesIO::writeFlowToFile()
{

}
