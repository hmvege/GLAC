#ifndef TEST_H
#define TEST_H

#include "tests/testsuite.h"
#include "performancetests.h"
#include "parallelization/communicator.h"

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
