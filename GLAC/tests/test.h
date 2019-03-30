/*!
 * \brief Header containing functions needed for running unit tests, integration tests and validation testing.
 *
 * \todo add validation testing using Jacks chroma config. Load it from python and compare in unit testing program.
 * \todo split into unit tests, integration tests and validation testing.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef TEST_H
#define TEST_H

#include "tests/testsuite.h"
#include "performancetests.h"
#include "parallelization/communicator.h"

#include "tests/unit_tests/actiontests.h"
#include "tests/unit_tests/complexoperations.h"
#include "tests/unit_tests/iotests.h"
#include "tests/unit_tests/latticeoperations.h"
#include "tests/unit_tests/observabletests.h"
#include "tests/unit_tests/randommatrixtests.h"
#include "tests/unit_tests/su2operations.h"
#include "tests/unit_tests/su3operations.h"
//#include "tests/unit_tests/testcore.h"

/*!
 * \brief runUnitTests function that runs unit tests. Exits once complete.
 * \param runTests if true, will run unit tests.
 */
void runUnitTests(bool runTests)
{
    if (runTests)
    {
//        TestSuite unitTester;
//        unitTester.runFullTestSuite();

        bool passed = true;
        if (Parallel::Communicator::getProcessRank() == 0) {
            for (int i = 0; i < 60; i++) cout << "=";
            cout << endl;
        }

        ActionTests actionTest = ActionTests();
        passed = passed & actionTest.runActionTests();

        ComplexOperations cmplxTest = ComplexOperations();
        passed = passed & cmplxTest.runComplexTests();

        IOTests IOTest = IOTests();
        passed = passed & IOTest.runIOTests();

        LatticeOperations latTest = LatticeOperations();
        passed = passed & latTest.runLatticeTests();

        ObservableTests obsTest = ObservableTests();
        passed = passed & obsTest.runObservableTests();

        RandomMatrixTests rndTest = RandomMatrixTests();
        passed = passed & rndTest.runRandomMatrixTests();

        SU2Operations SU2Test = SU2Operations();
        passed = passed & SU2Test.runSU2Tests();

        SU3Operations SU3Test = SU3Operations();
        passed = passed & SU3Test.runSU3Tests();

        if (Parallel::Communicator::getProcessRank() == 0) {
            for (int i = 0; i < 60; i++) cout << "=";
            cout << endl;
            if (passed) {
                cout << "SUCCESS: All tests passed." << endl;
            } else {
                cout << "FAILURE: One or more test(s) failed." << endl;
            }
            for (int i = 0; i < 60; i++) cout << "=";
            cout << endl;
        }

        Parallel::Communicator::setBarrier();
        Parallel::Communicator::freeMPIGroups();
        Parallel::Communicator::setBarrier();
        MPI_Finalize();
        exit(0);
    }
}

/*!
 * \brief runPerformanceTests will run performance tests on program. Exits once complete.
 * \param runTests if true, will run performance tests.
 */
void runPerformanceTests(bool runTests)
{
    if (runTests)
    {
        PerformanceTests performenceTester;
        performenceTester.run();
        Parallel::Communicator::setBarrier();
        Parallel::Communicator::freeMPIGroups();
        Parallel::Communicator::setBarrier();
        MPI_Finalize();
        exit(0);
    }
}


#endif // TEST_H
