#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#define CATCH_CONFIG_RUNNER  // This tells Catch that I will provide a main()
                             // function, and ensures Catch::Sessions API is
                             // available.
#include <mpi.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_case_info.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/reporters/catch_reporter_streaming_base.hpp>
#include <parallelization/communicator.h>

namespace
{
  // Define ANSI color codes
  const std::string RED = "\033[31m";
  const std::string GREEN = "\033[32m";
  const std::string YELLOW = "\033[33m";
  const std::string CYAN = "\033[36m";
  const std::string RESET = "\033[0m";
  const std::string BOLD = "\033[1m";
}  // namespace

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

    if (!assertionStats.assertionResult.isOk())
    {
      std::ostringstream oss;
      //  Prints failed process
      oss << "\n" << Catch::lineOfChars('*') << "\n";
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
      oss << Catch::lineOfChars('.') << "\n";
      oss << m_sectionStack.back().lineInfo << "\n";
      oss << Catch::lineOfChars('.') << "\n";

      // Prints failed test
      oss << assertionStats.assertionResult.getSourceInfo() << ": " << RED
          << "FAILED:" << RESET << "\n";
      oss << "  " << assertionStats.assertionResult.getExpression() << "\n";
      oss << "with expansion\n";
      oss << "  " << assertionStats.assertionResult.getExpandedExpression()
          << "\n";

      oss << "\n";
      m_failure_messages.push_back(oss.str());
    }
  }

  void testRunEnded(Catch::TestRunStats const& testRunStats) override
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Use MPI to determine if any process has failed
    int localFailedProcs = testRunStats.totals.testCases.failed > 0 ? 1 : 0;
    int globalFailedProcs = 0;
    MPI_Allreduce(&localFailedProcs, &globalFailedProcs, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    int localAssertions = testRunStats.totals.assertions.total();
    int globalAssertions = 0;
    MPI_Allreduce(&localAssertions, &globalAssertions, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    int localAssertionsPassed = testRunStats.totals.assertions.passed;
    int globalAssertionsPassed = 0;
    MPI_Allreduce(&localAssertionsPassed, &globalAssertionsPassed, 1, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);

    int localAssertionsFailed = testRunStats.totals.assertions.failed;
    int globalAssertionsFailed = 0;
    MPI_Allreduce(&localAssertionsFailed, &globalAssertionsFailed, 1, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);

    int localCases = testRunStats.totals.testCases.total();
    int globalCases = 0;
    MPI_Allreduce(&localCases, &globalCases, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    int localCasesPassed = testRunStats.totals.testCases.passed;
    int globalCasesPassed = 0;
    MPI_Allreduce(&localCasesPassed, &globalCasesPassed, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    int localCasesFailed = testRunStats.totals.testCases.failed;
    int globalCasesFailed = 0;
    MPI_Allreduce(&localCasesFailed, &globalCasesFailed, 1, MPI_INT, MPI_SUM,
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
      m_stream << Catch::lineOfChars('-') << "\n";

      // If there's a global failure, print a custom summary
      if (globalCasesFailed > 0)
      {
        m_stream << "There were failures in the test cases.\n";
        m_stream << Catch::lineOfChars('=') << "\n";
        m_stream << failure_message << "\n";
      }
      else
      {
        // Print the normal summary
        m_stream << "All tests passed.\n";
        m_stream << Catch::lineOfChars('=') << "\n\n";
      }

      m_stream << Catch::lineOfChars('=') << "\n";
      m_stream << "Total test cases run in parallel: "
               << testRunStats.totals.testCases.total() << "\n";

      int width_passed = std::to_string(globalAssertions).size() - 1;
      int width_failed = std::to_string(globalAssertionsFailed).size() - 1;

      // Printing test cases
      m_stream << "Test cases (total): " << std::setw(width_passed)
               << globalCases << " | "
               << m_colour->guardColour(Catch::Colour::Green)
               << globalCasesPassed
               << m_colour->guardColour(Catch::Colour::Green) << " passed"
               << " | " << std::setw(width_failed)
               << m_colour->guardColour(Catch::Colour::Red) << globalCasesFailed
               << m_colour->guardColour(Catch::Colour::Red) << " failed"
               << "\n";

      // Printing assertions
      m_stream << "Assertions (total): " << std::setw(width_passed)
               << globalAssertions << " | "
               << m_colour->guardColour(Catch::Colour::Green)
               << globalAssertionsPassed
               << m_colour->guardColour(Catch::Colour::Green) << " passed"
               << " | " << std::setw(width_failed)
               << m_colour->guardColour(Catch::Colour::Red)
               << globalAssertionsFailed
               << m_colour->guardColour(Catch::Colour::Red) << " failed"
               << "\n";

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
  Parallel::Communicator::init(argc, argv);

  // Prepare command line arguments for Catch
  std::vector<std::string> args(argv, argv + argc);
  const bool reporterSpecified = std::find(args.begin(), args.end(), "--reporter")
                                   != args.end()
                                 || std::find(args.begin(), args.end(), "-r")
                                      != args.end();
  if (!reporterSpecified)
  {
    args.push_back("--reporter");
    args.push_back(MPI_CATCH_REPORTER);
  }
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
    Parallel::Communicator::freeMPIGroups();
    MPI_Finalize();
    return returnCode;
  }

  int numFailed = session.run();

  // Turn local failures into a 0/1 flag
  int localFail = (numFailed > 0) ? 1 : 0;
  int globalFail = 0;
  MPI_Allreduce(&localFail, &globalFail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  Parallel::Communicator::freeMPIGroups();
  MPI_Finalize();

  return globalFail ? 1 : 0;
}
