#ifndef TEST_H
#define TEST_H

#include "tests/testsuite.h"
#include "performancetests.h"

void runUnitTests(bool runTests)
{
    if (runTests)
    {
        TestSuite unitTester;
        unitTester.runFullTestSuite();
        Parallel::Communicator::MPIExit("Unit tests complete.");
    }
}

void runPerformanceTests(bool runTests)
{
    if (runTests)
    {
        PerformanceTests performenceTester;
        performenceTester.run();
        Parallel::Communicator::MPIExit("Performance tester complete.");
    }
}


#endif // TEST_H
