#include "sysprint.h"
#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

SysPrint::SysPrint()
{
}

/*!
 * \brief SysPrint::printLine prints a simple line for seperating outputs.
 */
void SysPrint::printLine()
{
    for (int i = 0; i < 100; i++) cout << "=";
    cout << endl;
}

std::string SysPrint::getListString(const std::vector<std::string> &observableList)
{
    std::string returnString = observableList[0];
    for (unsigned int i = 1; i < observableList.size(); i++) {
        returnString += ", " + observableList[i];
    }
    return returnString;
}

std::string SysPrint::getTrueOrFalseString(const bool i)
{
    if (i) {
        return "TRUE";
    } else {
        return "FALSE";
    }
}

/*!
 * \brief SysPrint::printSystemInfo prints out parametesr from Parallel:Communicator and Parameters.
 */
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
        cout << "Sub lattice dimensions:                ";
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
        if (!Parameters::getLoadFieldConfigurations()) {
            cout << "N configurations:                      " << Parameters::m_NCf << endl;
            cout << "N correlation updates:                 " << Parameters::m_NCor << endl;
            cout << "N thermalization updates:              " << Parameters::m_NTherm << endl;
        }
        if (Parameters::m_NFlows != 0) {
            cout << "N flow updates per configuration:      " << Parameters::m_NFlows << endl;
        }
        if (!Parameters::getLoadFieldConfigurations()) {
            cout << "N link updates:                        " << Parameters::m_NUpdates << endl;
        }
        cout << "Output folder:                         " << Parameters::m_pwd + Parameters::m_outputFolder << endl;
        cout << "Input folder:                          " << Parameters::m_pwd + Parameters::m_inputFolder << endl;
        if (!Parameters::getLoadFieldConfigurations()) {
            cout << "Store field configurations:            " << getTrueOrFalseString(Parameters::m_storeConfigurations) << endl;
            cout << "Store thermalization observables:      " << getTrueOrFalseString(Parameters::m_storeThermalizationObservables) << endl;
        }
        if (Parameters::m_NFlows != 0) {
            cout << "Exponentiation function:               " << Parameters::m_expFuncName << endl;
        }
        cout << "Action type:                           " << Parameters::getActionType() << endl;
        if (!Parameters::getLoadFieldConfigurations()) {
            cout << "Observables:                           " << getListString(Parameters::m_observablesList) << endl;
        }
        if (Parameters::m_NFlows != 0) {
            cout << "Flow observables:                      " << getListString(Parameters::m_flowObservablesList) << endl;
            cout << "Flow epsilon:                          " << Parameters::m_flowEpsilon << endl;
            if (Parameters::getFlowObservablesList()[0] == "energyTopcFieldDensity")
            {
                cout << "Flow sampling frequency:               " << Parameters::getSamplingFrequency() << endl;
            }
        }
        if (!Parameters::getLoadFieldConfigurations()) {
            cout << "SU3Eps:                                " << Parameters::m_SU3Eps << endl;
            if (Parameters::m_verbose) {
                cout << "Metropolis seed:                       ";
                printf("%-16.1f\n",Parameters::m_metropolisSeed);
                cout << "Random matrix seed:                    ";
                printf("%-16.1f\n",Parameters::m_randomMatrixSeed);
            }
            cout << "Hot start:                             " << getTrueOrFalseString(Parameters::m_hotStart) << endl;
            if (Parameters::m_hotStart) {
                cout << "Random start close to unity:           " << getTrueOrFalseString(Parameters::m_RSTHotStart) << endl;
            }
        }
        if (Parameters::m_debug) {
            cout << "Debug:                                 TRUE" << endl;
        }
        printLine();
    }
}
