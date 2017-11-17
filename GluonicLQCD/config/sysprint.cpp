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

void SysPrint::printSystemInfo()
{
    /*
     * Function for printing system information in the beginning.
     */
    if (m_processRank == 0) {
        printLine();
        cout << "Batch name:                            " << Parameters::m_batchName << endl;
        cout << "Threads:                               " << Parallel::Communicator::m_numprocs << endl;
        cout << "Lattice size:                          " << Parameters::m_latticeSize << endl;
        cout << "Lattice dimensions(spatial, temporal): " << Parameters::m_NSpatial << " " << Parameters::m_NTemporal << endl;
        cout << "N configurations:                      " << Parameters::m_NCf << endl;
        cout << "N flow updates per configuration:      " << Parameters::m_NFlows << endl;
        cout << "N correlation updates:                 " << Parameters::m_NCor << endl;
        cout << "N thermalization updates:              " << Parameters::m_NTherm << endl;
        cout << "N link updates:                        " << Parameters::m_NUpdates << endl;
        cout << "Beta:                                  " << Parameters::m_beta << endl;
        cout << "SU3Eps:                                " << Parameters::m_SU3Eps << endl;
        cout << "Sub lattice Size:                      " << Parameters::m_subLatticeSize << endl;
        cout << "Sub latticedimensions:                 ";
        for (int i = 0; i < 4; i++) {
            cout << Parameters::m_N[i] << " ";
        }
        cout << endl;
        cout << "Processsors per dimension:             ";
        for (int i = 0; i < 4; i++) {
            cout << Parameters::m_processorsPerDimension[i] << " ";
        }
        cout << endl;
        printLine();
    }
}
