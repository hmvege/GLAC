#include "observabletests.h"
#include "observables/observables.h"
#include "parallelization/communicator.h"

ObservableTests::ObservableTests()
{

}


bool ObservableTests::testTopCharge()
{
    bool passed = true;

    return passed;
}

// Observable tests
bool ObservableTests::runObservableTests()
{
    bool passed = false;

    if (m_processRank == 0) {
        if (m_verbose) cout << "Running observables tests." << endl;

        passed = testTopCharge();

        if (passed) {
            cout << "PASSED: observable tests." << endl;
        } else {
            cout << "FAILURE: observable tests." << endl;
        }
    }

    MPI_Bcast(&passed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    Parallel::Communicator::setBarrier();

    return passed;
}
