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
  std::vector<std::string> m_failure_messages;

public:
  using StreamingReporterBase::StreamingReporterBase;

  void assertionEnded(Catch::AssertionStats const& assertionStats) override
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::cout << "Rank: " << rank << "\n";

    if (!assertionStats.assertionResult.isOk())
    {
      // m_stream << StreamingReporterBase::currentTestCaseInfo->name << "\n";
      // m_stream << StreamingReporterBase::currentTestCaseInfo->className <<
      // "\n"; m_stream << StreamingReporterBase::currentTestRunInfo.name <<
      // "\n";

      // m_stream << Catch::lineOfChars('-') << "\n";
      // m_stream << "Assertion failed at process of: "
      //          << m_colour->guardColour(Catch::Colour::Yellow) << "RANK "
      //          << rank << "\n";

      // // Printing source file and line
      // m_stream << assertionStats.assertionResult.getSourceInfo() << ": ";
      // m_stream << m_colour->guardColour(Catch::Colour::Red) << "FAILED";
      // // Printing out the expression that failed.
      // m_stream << "\n\t"
      //          << assertionStats.assertionResult.getExpressionInMacro() <<
      //          "\n";

      // // In case there are any additional messages(e.g. captures), we print
      // // those.
      // if (assertionStats.assertionResult.hasMessage())
      // {
      //   m_stream << assertionStats.assertionResult.getMessage() << "\n";
      // }
      // m_stream << Catch::lineOfChars('.') << "\n";

      std::ostringstream oss;
      //  Prints failed process
      oss << Catch::lineOfChars('*') << "\n";
      oss << "Assertion failed at process RANK " << rank << ":\n";

      // Prints location of the failed test
      oss << Catch::lineOfChars('-') << "\n";
      if (m_sectionStack.size() > 0)
      {
        oss << m_sectionStack.at(0).name << "\n";
      }
      if (m_sectionStack.size() > 1)
      {
        for (int i = 1; i < m_sectionStack.size(); i++)
        {
          oss << "  " << m_sectionStack.at(i).name << "\n";
        }
      }

      // Prints section location
      oss << Catch::lineOfChars('-') << "\n";
      oss << m_sectionStack.back().lineInfo << "\n";
      oss << Catch::lineOfChars('.') << "\n";

      // Prints failed test
      oss << assertionStats.assertionResult.getSourceInfo() << ": FAILED:\n";
      oss << "  " << assertionStats.assertionResult.getExpression() << "\n";
      oss << "with expansion\n";
      oss << "  " << assertionStats.assertionResult.getExpandedExpression()
          << "\n";

      // oss << "  Test case: " << currentTestCaseInfo->name << "\n";
      // oss << "  Assertion: "
      //     << assertionStats.assertionResult.getExpressionInMacro() << "\n";
      // oss << "  Message: " << assertionStats.assertionResult.getMessage()
      //     << "\n";

      // m_stream << "getExpandedExpression: "
      //          << assertionStats.assertionResult.getExpandedExpression()
      //          << "\n";  // Fetches full expression! USE THIS!
      // m_stream << "  Message: " <<
      // assertionStats.assertionResult.getMessage()
      //          << "\n";
      // m_stream << "getExpressionInMacro: "
      //          << assertionStats.assertionResult.getExpressionInMacro() <<
      //          "\n";
      // m_stream << "getExpression: "
      //          << assertionStats.assertionResult.getExpression() << "\n";
      // m_stream << "getSourceInfo: "
      //          << assertionStats.assertionResult.getSourceInfo() << "\n";
      // m_stream << "getTestMacroName: "
      //          << assertionStats.assertionResult.getTestMacroName() << "\n";

      oss << Catch::lineOfChars('-') << "\n";
      m_failure_messages.push_back(oss.str());
    }
    // StreamingReporterBase::assertionEnded(assertionStats);

    m_stream << "ASSERTION END\n";
  }

  // // Override necessary methods
  // void testCaseEnded(Catch::TestCaseStats const& testCaseStats) override
  // {
  //   int rank;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //   // MPI_Barrier(MPI_COMM_WORLD);

  //   if (testCaseStats.totals.testCases.failed > 0)
  //   {
  //     m_stream << "\n";
  //     m_stream << "Test case failed at process of: "
  //              << m_colour->guardColour(Catch::Colour::Yellow) << "RANK "
  //              << rank << "\n";
  //     m_stream << "Failed test: " << testCaseStats.testInfo->name << "\n";
  //     // testCaseStats.testInfo->properties
  //     // for (const auto t : testCaseStats.testInfo->tags)
  //     // {
  //     //   m_stream << t.original << "\n";
  //     // }
  //     m_stream << "    Failed at: " << testCaseStats.testInfo->lineInfo <<
  //     "\n"; if (testCaseStats.stdOut.length() > 0)
  //     {
  //       m_stream << testCaseStats.stdOut << "\n";
  //     }
  //     if (testCaseStats.stdErr.length() > 0)
  //     {
  //       m_stream << testCaseStats.stdErr << "\n";
  //     }
  //     // for (auto c : testCaseStats.testInfo->tags)
  //     //   std::cout << c.original
  //     // m_stream << testCaseStats.totals.assertions.failed << "\n";
  //     m_stream << Catch::lineOfChars('-') << "\n";
  //   }
  //   StreamingReporterBase::testCaseEnded(testCaseStats);
  // }

  // void sectionEnded(Catch::SectionStats const& sectionStats)
  // {
  //   int rank;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   if (!sectionStats.assertions.allPassed())
  //   {
  //     // m_stream << "Test section failed at process of: "
  //     //          << m_colour->guardColour(Catch::Colour::Yellow) << "RANK "
  //     //          << rank << "\n";
  //     m_stream << "\t" << sectionStats.sectionInfo.name << "\n";
  //     m_stream << Catch::lineOfChars('#') << "\n";
  //   }
  //   StreamingReporterBase::sectionEnded(sectionStats);
  // }

  void testRunEnded(Catch::TestRunStats const& testRunStats) override
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Use MPI to determine if any process has failed
    int localFailedProcs = testRunStats.totals.testCases.failed > 0 ? 1 : 0;
    int globalFailedProcs = 0;
    MPI_Allreduce(&localFailedProcs, &globalFailedProcs, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    int totalLocalPassed = testRunStats.totals.testCases.passed;
    int totalGlobalPassed = 0;
    MPI_Allreduce(&totalLocalPassed, &totalGlobalPassed, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    int totalLocalFailed = testRunStats.totals.testCases.failed;
    int totalGlobalFailed = 0;
    MPI_Allreduce(&totalLocalFailed, &totalGlobalFailed, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    int totalLocalSkipped = testRunStats.totals.testCases.skipped;
    int totalGlobalSkipped = 0;
    MPI_Allreduce(&totalLocalSkipped, &totalGlobalSkipped, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    // Fetches all of the local messages into 1 long string
    std::string localFailuresMessage;
    for (const auto msg : m_failure_messages)
    {
      localFailuresMessage += msg;
    }
    // Fetches size of all the local messages
    const int localMessageSize = localFailuresMessage.size();

    // Gather the size of all messages for rank 0 processor
    std::vector<int> recvSizes(size);
    MPI_Gather(&localMessageSize, 1, MPI_INT, recvSizes.data(), 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    // Now, we can gather all of the failure messages into a single long string.
    // First, we calculate the total size for rank 0 and size
    int totalGlobalSize = 0;
    std::vector<int> recvDisplacements(size);
    if (rank == 0)
    {
      for (int i = 0; i < recvSizes.size(); i++)
      {
        recvDisplacements.at(i) = totalGlobalSize;
        totalGlobalSize += recvSizes.at(i);
      }
    }

    // Allocate space for messages
    std::vector<char> allFailures(totalGlobalSize);
    MPI_Gatherv(localFailuresMessage.data(), localMessageSize, MPI_CHAR,
                allFailures.data(), recvSizes.data(), recvDisplacements.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD);

    std::string failure_message(allFailures.begin(), allFailures.end());

    // Only the master process prints the summary
    if (rank == 0)
    {
      m_stream << Catch::lineOfChars('=') << "\n";

      m_stream << "Catch tests for '" << testRunStats.runInfo.name
               << "' completed.\n";
      m_stream << Catch::lineOfChars('=') << "\n";

      // If there's a global failure, print a custom summary
      if (totalGlobalFailed > 0)
      {
        m_stream << "There were failures in the test cases.\n";
        m_stream << failure_message << "\n";
      }
      else
      {
        // Print the normal summary
        m_stream << "All tests passed.\n";
      }

      m_stream << Catch::lineOfChars('=') << "\n";
      m_stream << "Total test cases: " << testRunStats.totals.testCases.total()
               << "\n";
      m_stream << m_colour->guardColour(Catch::Colour::Green)
               << "PASSED: " << totalGlobalPassed << "\n";
      m_stream << m_colour->guardColour(Catch::Colour::Red)
               << "FAILED: " << totalGlobalFailed << "\n";

      if (totalGlobalSkipped > 0)
      {
        m_stream << m_colour->guardColour(Catch::Colour::Skip)
                 << "SKIPPED: " << totalGlobalSkipped << "\n";
      }

      m_stream << Catch::lineOfChars('=') << "\n";
    }
    StreamingReporterBase::testRunEnded(testRunStats);
  }

  static std::string getDescription() { return "MPI Custom Catch reporter"; }
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