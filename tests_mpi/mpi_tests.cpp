#include <lattice/test_lattice.h>
#include <parallelization/communicator.h>

#include <iomanip>
#include <iostream>

int main(int numberOfArguments, char* cmdLineArguments[])
{
  Parallel::Communicator::init(numberOfArguments, cmdLineArguments);

  std::cout << "Starting MPI tests\n";

  bool passed = true;

  if (Parallel::ParallelParameters::active)
  {
    // Test math
    passed &= testLattice();

    // TODO: implement this properly after clean-up of observables has been
    // performed.
    // Test observables
    // passed &= testCorrelator();
    // passed &= testObservablesStorer();

    // System actions
    // passed &= testAction();
    // passed &= testWilsonExplicitDer();
    // passed &= testGaugeAction();

    // System flow
    // passed &= testFlow();

    // System tests
    // passed &= testSystem();

    // Config tests
    // passed &= testParameters(); // TODO: implement after improved parameters
    // passed &= testConfigLoader();
    // passed &= testSystemPrint();

    // Parallelization tests
    // passed &= testIndex(); // TODO: perhaps use MPIs implementation
    // passed &= testParallelParameters();
    // passed &= testCommunicator();
    // passed &= testNeighbourList();
    // passed &= testNeighbours();
  }

  return passed && MPI_Finalize();
}

// #define CATCH_CONFIG_RUNNER  // This tells Catch to provide a main() function
// #include <mpi.h>

// #include <catch2/catch_all.hpp>

// int main(int argc, char* argv[])
// {
//   MPI_Init(&argc, &argv);

//   Catch::Session session;

//   int returnCode = session.applyCommandLine(argc, argv);

//   if (returnCode != 0)
//   {
//     MPI_Finalize();
//     return returnCode;
//   }

//   int numFailed = session.run();

//   MPI_Finalize();

//   return numFailed;
// }

// class PartialReporter : public Catch::StreamingReporterBase {
// public:
// void ConsoleReporter::testRunEnded(TestRunStats const& _testRunStats)
// {
//   int rank id = -1;
//   MPI Comm rank(MPI COMM WORLD, &rank id);
//   if (rank id != 0 && testRunStats.totals.testCases.allPassed()) return;
//   printTotalsDivider(_testRunStats.totals);
//   printTotals(_testRunStats.totals);
//   stream << std::endl;
//   StreamingReporterBase::testRunEnded(_testRunStats);
// }
// }

// #include <catch2/catch_test_case_info.hpp>
// #include <catch2/reporters/catch_reporter_registrars.hpp>
// #include <catch2/reporters/catch_reporter_streaming_base.hpp>
// #include <iostream>

// class PartialReporter : public Catch::StreamingReporterBase
// {
// public:
//   using StreamingReporterBase::StreamingReporterBase;

//   static std::string getDescription()
//   {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0)
//       return "Reporter for testing TestCasePartialStarting/Ended events";
//     else
//       return "";
//   }

//   void testCasePartialStarting(Catch::TestCaseInfo const& testInfo,
//                                uint64_t partNumber) override
//   {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0)
//       std::cout << "TestCaseStartingPartial: " << testInfo.name << '#'
//                 << partNumber << '\n';
//   }

//   void testCasePartialEnded(Catch::TestCaseStats const& testCaseStats,
//                             uint64_t partNumber) override
//   {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (rank == 0)
//       std::cout << "TestCasePartialEnded: " << testCaseStats.testInfo->name
//                 << '#' << partNumber << '\n';
//   }
// };

// // CATCH_REGISTER_REPORTER("partial", PartialReporter)