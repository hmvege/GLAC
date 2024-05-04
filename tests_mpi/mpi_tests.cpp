// #include <lattice/test_lattice.h>
// #include <parallelization/communicator.h>

// #include <iomanip>
// #include <iostream>

// int main(int numberOfArguments, char* cmdLineArguments[])
// {
//   Parallel::Communicator::init(numberOfArguments, cmdLineArguments);

//   std::cout << "Starting MPI tests\n";

//   bool passed = true;

//   if (Parallel::ParallelParameters::active)
//   {
//     // Test math
//     // passed &= testLattice();

//     // TODO: implement this properly after clean-up of observables has been
//     // performed.
//     // Test observables
//     // passed &= testCorrelator();
//     // passed &= testObservablesStorer();

//     // System actions
//     // passed &= testAction();
//     // passed &= testWilsonExplicitDer();
//     // passed &= testGaugeAction();

//     // System flow
//     // passed &= testFlow();

//     // System tests
//     // passed &= testSystem();

//     // Config tests
//     // passed &= testParameters(); // TODO: implement after improved
//     parameters
//     // passed &= testConfigLoader();
//     // passed &= testSystemPrint();

//     // Parallelization tests
//     // passed &= testIndex(); // TODO: perhaps use MPIs implementation
//     // passed &= testParallelParameters();
//     // passed &= testCommunicator();
//     // passed &= testNeighbourList();
//     // passed &= testNeighbors();
//   }

//   return passed && MPI_Finalize();
// }

#define CATCH_CONFIG_RUNNER  // This tells Catch to provide a main() function
#include <mpi.h>

#include <catch2/catch_all.hpp>
// #include <catch2/catch_reporter_bases.hpp>
#include <catch2/catch_test_case_info.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/reporters/catch_reporter_streaming_base.hpp>

namespace
{}  // namespace

class MPIMasterReporter : public Catch::StreamingReporterBase
{
private:
public:
  using StreamingReporterBase::StreamingReporterBase;

  void assertionEnded(Catch::AssertionStats const& assertionStats)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!assertionStats.assertionResult.isOk())
    {
      m_stream << StreamingReporterBase::currentTestCaseInfo->name << "\n";
      m_stream << StreamingReporterBase::currentTestCaseInfo->className << "\n";
      m_stream << StreamingReporterBase::currentTestRunInfo.name << "\n";

      m_stream << Catch::lineOfChars('-') << "\n";
      m_stream << "Assertion failed at process of: "
               << m_colour->guardColour(Catch::Colour::Yellow) << "RANK "
               << rank << "\n";

      // Printing source file and line
      m_stream << assertionStats.assertionResult.getSourceInfo() << ": ";
      m_stream << m_colour->guardColour(Catch::Colour::Red) << "FAILED";
      // Printing out the expression that failed.
      m_stream << "\n\t"
               << assertionStats.assertionResult.getExpressionInMacro() << "\n";

      // In case there are any additional messages(e.g. captures), we print
      // those.
      if (assertionStats.assertionResult.hasMessage())
      {
        m_stream << assertionStats.assertionResult.getMessage() << "\n";
      }

      m_stream << Catch::lineOfChars('.') << "\n";
    }
    StreamingReporterBase::assertionEnded(assertionStats);

    // m_stream << "ASSERTION END\n";
  }

  // Override necessary methods
  void testCaseEnded(Catch::TestCaseStats const& testCaseStats) override
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if (testCaseStats.totals.testCases.failed > 0)
    {
      m_stream << "Test case failed at process of: "
               << m_colour->guardColour(Catch::Colour::Yellow) << "RANK "
               << rank << "\n";
      m_stream << testCaseStats.stdOut << "\n";
      m_stream << testCaseStats.stdErr << "\n";
      m_stream << testCaseStats.testInfo->lineInfo << "\n";
      m_stream << testCaseStats.testInfo->name << "\n";
      // for (auto c : testCaseStats.testInfo->tags)
      //   std::cout << c.original
      for (const auto t : testCaseStats.testInfo->tags)
        m_stream << t.original << "\n";
      m_stream << testCaseStats.totals.assertions.failed << "\n";
      m_stream << Catch::lineOfChars('%') << "\n";
    }
    StreamingReporterBase::testCaseEnded(testCaseStats);
  }

  void sectionEnded(Catch::SectionStats const& sectionStats)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if (!sectionStats.assertions.allPassed())
    {
      // m_stream << "Test section failed at process of: "
      //          << m_colour->guardColour(Catch::Colour::Yellow) << "RANK "
      //          << rank << "\n";
      m_stream << "\t" << sectionStats.sectionInfo.name << "\n";
      m_stream << Catch::lineOfChars('#') << "\n";
    }
    StreamingReporterBase::sectionEnded(sectionStats);
  }

  void testRunEnded(Catch::TestRunStats const& testRunStats) override
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Use MPI to determine if any process has failed
    int localFail = testRunStats.totals.testCases.failed > 0 ? 1 : 0;
    int globalFail;
    MPI_Allreduce(&localFail, &globalFail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Waits for all processes to finish, then print results.
    MPI_Barrier(MPI_COMM_WORLD);

    // Only the master process prints the summary
    if (rank == 0)
    {
      m_stream << Catch::lineOfChars('=') << "\n";

      m_stream << "Catch tests for '" << testRunStats.runInfo.name
               << "' completed.\n";
      m_stream << Catch::lineOfChars('.') << "\n";

      m_stream << "Total test cases: " << testRunStats.totals.testCases.total()
               << "\n";
      m_stream << m_colour->guardColour(Catch::Colour::Green)
               << "PASSED: " << testRunStats.totals.assertions.passed << "\n";
      m_stream << m_colour->guardColour(Catch::Colour::Red)
               << "FAILED: " << testRunStats.totals.assertions.failed << "\n";

      if (const auto skipped = testRunStats.totals.assertions.skipped;
          skipped > 0)
      {
        m_stream << m_colour->guardColour(Catch::Colour::Skip)
                 << "SKIPPED: " << testRunStats.totals.assertions.skipped
                 << "\n";
      }

      // If there's a global failure, print a custom summary
      if (globalFail > 0)
      {
        m_stream << "There were failures in the test cases.\n";
      }
      else
      {
        // Print the normal summary
        m_stream << "All tests passed.\n";
      }
      m_stream << Catch::lineOfChars('.') << "\n";
    }
    StreamingReporterBase::testRunEnded(testRunStats);
  }

  static std::string getDescription() { return {"MPI Custom Catch reporter"}; }
};

const std::string MPI_CATCH_REPORTER{"mpi_catch_reporter"};

CATCH_REGISTER_REPORTER(MPI_CATCH_REPORTER, MPIMasterReporter)

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  // Prepare command line arguments for Catch
  std::vector<std::string> args(argv, argv + argc);
  args.push_back("--reporter");
  args.push_back(MPI_CATCH_REPORTER);
  args.push_back("--colour-mode");
  args.push_back("ansi");

  // Convert modified arguments back to char* array
  std::vector<char*> argv_modified;
  for (auto& arg : args)
  {
    argv_modified.push_back(&arg[0]);
  }

  Catch::Session session;

  int returnCode = session.applyCommandLine(
    static_cast<int>(argv_modified.size()), argv_modified.data());

  if (returnCode != 0)
  {
    MPI_Finalize();
    return returnCode;
  }

  // TODO: make a dummy test which tests with proc 1 and 3 failing ,and see that
  // printing if happening in order

  int numFailed = session.run();

  MPI_Finalize();

  return 0;
}