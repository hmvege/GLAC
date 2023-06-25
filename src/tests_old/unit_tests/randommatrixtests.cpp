#include "randommatrixtests.h"

#include "parallelization/communicator.h"

RandomMatrixTests::RandomMatrixTests() {}

bool RandomMatrixTests::testRSTMultiplication()
{
  /*
   * Function for checking that we are multiplying the SU2 matrices correctly
   * when generating a SU3 matrix close to unity.
   */
  bool passed = true;
  SU3 results = m_SU3Generator->RSTMatrixMultiplication(s_r, s_s, s_t);
  if (compareSU3(results, U_RST))
  {
    if (m_verbose && (m_processRank == 0))
      cout << "    SUCCESS: RST multiplication test passed." << endl;
  }
  else
  {
    if (m_processRank == 0)
    {
      if (m_verbose) results.print();
      cout << "    FAILED: RST multiplication test did not pass." << endl;
    }
    passed = false;
  }
  return passed;
}

bool RandomMatrixTests::testRSTInverseMultiplication()
{
  /*
   * Function for checking that we are multiplying the SU2 matrices correctly
   * when generating a SU3 matrix close to unity.
   */
  bool passed = true;
  SU3 results = m_SU3Generator->RSTMatrixMultiplicationInverse(s_r, s_s, s_t);
  if (compareSU3(results, U_RST.inv()))
  {
    if (m_verbose)
      cout << "    SUCCESS: RST inverse multiplication test passed." << endl;
  }
  else
  {
    if (m_verbose) results.print();
    cout << "    FAILED: RST inverse multiplication test did not pass." << endl;
    passed = false;
  }
  return passed;
}

bool RandomMatrixTests::runRandomMatrixTests()
{
  bool passed = false;
  if (m_processRank == 0)
  {
    if (m_verbose)
      cout << "Running tests for Random matrix generation." << endl;

    passed = testRSTMultiplication() && testRSTInverseMultiplication();

    if (passed)
    {
      cout << "PASSED: random matrix generation methods." << endl;
    }
    else
    {
      cout << "FAILURE: random matrix generation methods." << endl;
    }
  }

  MPI_Bcast(&passed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
  Parallel::Communicator::setBarrier();

  return passed;
}
