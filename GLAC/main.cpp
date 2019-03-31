#include <chrono>
#include "system.h"
#include "config/parameters.h"
#include "config/configloader.h"
#include "tests/test.h"

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

/*!
 * \mainpage Documentation for GLAC
 *
 * GLAC is a program which generates gauge configurations and has the possibility of applying gradient flow on them as well.
 *
 * It takes input in the form of .json files generated by createJobs.py, a python2 command line interfacing tool. createJobs.py will create relevant folders for the run.
 * In the case that we are running on an HPC cluster using Slurm or Torque, it will directly initiate the run by submitting the jobs unless <code>--dryrun</code> is specified.
 *
 * See \ref page1 "Quickstart to GLAC" for a short guide in how to use GLAC.
 *
 * \sa
 * - <a href="https://github.com/hmvege/LQCDMasterThesis">The M.Sc. thesis by Mathias M Vege which was the basis for this code.</a>
 * - <a href="https://github.com/hmvege/LatViz">LatViz: A lattice visulization program</a> which can generate animations based on the output of the LatticeActionChargeDensity observable.
 */

/*!
 * \page page1 Quickstart to GLAC.
 * \code{.unparsed}
 *      python2 createJobs.py setup local 8 -N 16 -NT 32 -fobs topct -NFlows 1000 -NCfgs 100 -NTh 100000 -fEps 0.01 -b 6.0 -rn GLACTestRun
 * \endcode
 * This will generate a .json file called <code>config_GLACTestRun.json</code> located in a folder called <code>input</code>. The system will be of size \f$ 16^3 \times 32\f$,
 * and have \f$\beta = 6.0\f$. It will run \f$N_\mathrm{th} = 10^5\f$ thermalization steps, generate \f$N_\mathrm{cfg} = 100\f$ configurations, and flow each of the $\f$N_\mathrm{flow} = 1000\f$ a
 * steps with a step size of \f$\epsilon_\mathrm{flow} = 0.01\f$. The name of the run will be GLACTestRun.
 *
 * Type
 * \code{.unparsed}
 *      python2 createJobs.py -h
 * \endcode
 * in order to see what options are available.
 */


/*!
 * \brief main loads in configurations and initiates program according to .json specifications.
 *
 * \param numberOfArguments
 * \param cmdLineArguments
 * \return 0
 */
int main(int numberOfArguments, char* cmdLineArguments[])
{
    Parallel::Communicator::init(&numberOfArguments, &cmdLineArguments);

    ConfigLoader::load(std::string(cmdLineArguments[1]));

    // Initialises lattice sharing in the communicator
//    Parallel::Communicator::initializeSubLattice();

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

    // Finalizing and printing time taken
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\nProgram complete. Time used: %f hours (%f seconds)", double(programTime.count())/3600.0, programTime.count());
    }

    Parallel::Communicator::freeMPIGroups();
    Parallel::Communicator::setBarrier();
    MPI_Finalize();

    return 0;
}
