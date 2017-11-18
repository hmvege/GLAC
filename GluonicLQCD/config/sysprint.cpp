#include "sysprint.h"

using std::cout;
using std::endl;

SysPrint::SysPrint()
{
}

void SysPrint::printLine()
{
    for (int i = 0; i < 60; i++) cout << "=";
    cout << endl;
}

std::string SysPrint::getListString(std::vector<std::string> observableList)
{
    std::string returnString;
    for (unsigned int i = 0; i < observableList.size(); i++) {
        returnString += observableList[i];
    }
    return returnString;
}

std::string SysPrint::getTrueOrFalseString(bool i)
{
    if (i) {
        return "TRUE";
    } else {
        return "FALSE";
    }
}

void SysPrint::printSystemInfo()
{
    /*
     * Function for printing system information in the beginning.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        printLine();
        cout << "Batch name:                            " << Parameters::m_batchName << endl;
        cout << "Threads:                               " << Parallel::Communicator::getNumProc() << endl;
        cout << "Lattice size:                          " << Parameters::m_latticeSize << endl;
        cout << "Lattice dimensions(spatial, temporal): " << Parameters::m_NSpatial << " " << Parameters::m_NTemporal << endl;
        cout << "Sub lattice size:                      " << Parameters::m_subLatticeSize << endl;
        cout << "Sub lattice dimensions:                    ";
        for (int i = 0; i < 4; i++) {
            cout << Parameters::m_N[i] << " ";
        }
        cout << endl;
        cout << "Processsors per dimension:             ";
        for (int i = 0; i < 4; i++) {
            cout << Parameters::m_processorsPerDimension[i] << " ";
        }
        cout << endl;
        cout << "Beta:                                  " << Parameters::m_beta << endl;
        cout << "N configurations:                      " << Parameters::m_NCf << endl;
        cout << "N correlation updates:                 " << Parameters::m_NCor << endl;
        cout << "N thermalization updates:              " << Parameters::m_NTherm << endl;
        if (Parameters::m_NFlows != 0) {
            cout << "N flow updates per configuration:      " << Parameters::m_NFlows << endl;
        }
        cout << "N link updates:                        " << Parameters::m_NUpdates << endl;
        cout << "Output folder:                         " << Parameters::m_pwd + "/" + Parameters::m_outputFolder << endl;
        cout << "Input folder:                          " << Parameters::m_pwd + "/" + Parameters::m_inputFolder << endl;
        cout << "Store configurations:                  " << getTrueOrFalseString(Parameters::m_storeConfigurations) << endl;
        cout << "Store thermalization observables:      " << getTrueOrFalseString(Parameters::m_storeThermalizationObservables) << endl;
        cout << "Exponentiation function:               " << Parameters::m_expFuncName << endl;
        cout << "Observables:                           " << getListString(Parameters::m_observablesList) << endl;
        if (Parameters::m_NFlows != 0) {
            cout << "Flow observables:                      " << getListString(Parameters::m_flowObservablesList) << endl;
        }
        if (Parameters::m_NFlows != 0) {
            cout << "Flow epsilon:                          " << Parameters::m_flowEpsilon << endl;
        }
        cout << "SU3Eps:                                " << Parameters::m_SU3Eps << endl;
        if (Parameters::m_verbose) {
            cout << "Metropolis seed:                       " << Parameters::m_metropolisSeed << endl;
            cout << "Random matrix seed:                    " << Parameters::m_randomMatrixSeed << endl;
        }
        cout << "Hot start:                                 " << getTrueOrFalseString(Parameters::m_hotStart) << endl;

        printLine();
    }
}
