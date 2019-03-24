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

/*!
 * \brief runUnitTests function that runs unit tests. Exits once complete.
 * \param runTests if true, will run unit tests.
 */
void runUnitTests(bool runTests)
{
    if (runTests)
    {
        TestSuite unitTester;
        unitTester.runFullTestSuite();
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
