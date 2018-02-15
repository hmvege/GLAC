#include <chrono>
#include "system.h"
#include "config/parameters.h"
#include "config/configloader.h"
#include "tests/test.h"

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

void runUnitTests();
void runPerformanceTests();

int main(int numberOfArguments, char* cmdLineArguments[])
{
    Parallel::Communicator::init(&numberOfArguments, &cmdLineArguments);
    ConfigLoader::load(std::string(cmdLineArguments[1]));

    // Unit and performance tests
    runUnitTests(Parameters::getUnitTesting());
    runPerformanceTests(Parameters::getPerformanceTesting());

    // Program timers
    steady_clock::time_point programStart;
    programStart = steady_clock::now();

    // Main program part
    System pureGauge;
    pureGauge.latticeSetup();
    pureGauge.run();

//    Parallel::Communicator::MPIPrint("OK AFTER PURE GAUGE RUN?");

    // Finalizing and printing time taken
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    char msg[150];
    sprintf(msg,"\nProgram complete. Time used: %f hours (%f seconds)", double(programTime.count())/3600.0, programTime.count());
    Parallel::Communicator::MPIExit(std::string(msg));

    return 0;
}
