#include "su2operations.h"

#include <iomanip>

#include "parallelization/communicator.h"

SU2Operations::SU2Operations() {}

// Matrix operation tester
bool SU2Operations::operationSU2Test(SU2 results, SU2 solution,
                                     std::string operation)
{
  /*
   * Function that checks if a opertion on a 2x2 matrix has succeded.
   */
  if (compareSU2(results, solution))
  {
    if (m_verbose)
    {
      cout << "    SUCCESS: Matrix " << operation << " passed" << endl;
    }
    return true;
  }
  else
  {
    cout << "    FAILED: Matrix " << operation << " failed" << endl;
    if (m_verbose)
    {
      cout << "Results = " << endl;
      results.print();
      cout << "Solution = " << endl;
      solution.print();
    }
    return false;
  }
}

////////////////////////////////////
///////// 2x2 MATRIX TESTS /////////
////////////////////////////////////
bool SU2Operations::testSU2Addition()
{
  return operationSU2Test(s1 + s2, sAdd, "addition");
}

bool SU2Operations::testSU2Subtraction()
{
  return operationSU2Test(s1 - s2, sSub, "subtraction");
}

bool SU2Operations::testSU2Multiplication()
{
  return operationSU2Test(s1 * s2, sMul, "multiplication");
}

bool SU2Operations::testSU2Transpose()
{
  s3 = s1;
  return operationSU2Test(s3.transpose(), sT, "transpose");
}

bool SU2Operations::testSU2Conjugation()
{
  s3 = s1;
  return operationSU2Test(s3.conjugate(), sC, "conjugate");
}

bool SU2Operations::testSU2ComplexConjugation()
{
  s3 = s1;
  s3.conjugate();
  return operationSU2Test(s3.transpose(), sCT, "conjugate transpose");
}

////////////////////////////////////
/////// SU2 PROPERTIES TESTS ///////
////////////////////////////////////
bool SU2Operations::testSU2Hermicity()
{
  /*
   * Overarching function for checking the hermicity of the SU2 matrix.
   */
  bool passed = true;
  if (!checkSU2Hermicity(m_SU3Generator->generateSU2()))
  {
    cout << "    FAILED: randomly generated SU2 matrix is not hermitian."
         << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: randomly generated SU2 matrix is hermitian." << endl;
  }
  return passed;
}

bool SU2Operations::checkSU2Hermicity(SU2 H)
{
  /*
   * Checks the hermicity of the SU2 matrices.
   * H = 0 1  2 3
   *     4 5  6 7
   */
  bool passed = true;
  SU2 I;
  I = H * H.inv();
  for (int i = 0; i < 8; i++)
  {
    if (i == 0 || i == 6) continue;
    if (fabs(I[i]) > m_eps) passed = false;
  }
  if ((fabs(I[0] - 1) > m_eps) || (fabs(I[6] - 1) > m_eps)) passed = false;

  if (passed)
  {
    return true;
  }
  else
  {
    if (m_verbose)
    {
      cout << "H = " << endl;
      H.print();
      cout << "H^-1 = " << endl;
      H.inv().print();
      cout << "H^-1*H = " << endl;
      I.print();
    }
    return false;
  }
}

bool SU2Operations::testSU2Orthogonality()
{
  /*
   * Checks if SU2 matrices generated randomly by RST or Random is orthogonal.
   */
  bool passed = true;
  if (!checkSU2Orthogonality(m_SU3Generator->generateSU2()))
  {
    cout << "    FAILED: SU2 randomly generated matrix is not orthogonal."
         << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: randomly generated SU2 matrix is orthogonal." << endl;
  }
  return passed;
}

bool SU2Operations::checkSU2Orthogonality(SU2 H)
{
  /*
   * Testing the orthogonatility of a SU3 matrix H.
   */
  bool passed = true;
  complex col1[2];
  complex col2[2];
  for (int i = 0; i < 2; i++)
  {
    col1[i].setRe(H[4 * i]);
    col1[i].setIm(H[4 * i + 1]);
    col2[i].setRe(H[4 * i + 2]);
    col2[i].setIm(H[4 * i + 3]);
  }
  complex c12dot = dot2(col1, col2);
  if ((fabs(c12dot.re()) > m_eps) || (fabs(c12dot.im()) > m_eps))
    passed = false;
  if (m_verbose && !passed)
  {
    cout << "Column 1 and 2: " << std::setprecision(10) << c12dot << endl;
  }
  return passed;
}

bool SU2Operations::testSU2Norm()
{
  /*
   * Overarching function for testing if the columns of the SU3 matrices
   * generated by the RST and Random method are at unity.
   */
  bool passed = true;
  if (!checkSU2Norm(m_SU3Generator->generateSU2()))
  {
    cout << "    FAILED: columns of a randomly generated SU2 matrix is not at "
            "length 1."
         << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: randomly generated SU2 matrices have columns of "
            "length 1."
         << endl;
  }
  return passed;
}

bool SU2Operations::checkSU2Norm(SU2 H)
{
  /*
   * Function that checks the SU3 norm of a matrix H.
   */
  bool passed = true;
  double sum;
  for (int col = 0; col < 2; col++)
  {
    sum = 0;
    for (int i = 0; i < 2; i++)
    {
      sum += H.normSquared(4 * i + 2 * col);
    }
    sum = sqrt(sum);
    if (fabs(sqrt(sum) - 1) > m_eps)
    {
      if (m_verbose) cout << "    FAILED: norm is not zero: " << sum << endl;
      return false;
    }
  }
  return passed;
}

bool SU2Operations::testSU2Determinant()
{
  /*
   * Overarching function for testing the SU2 determinant.
   */
  bool passed = true;
  if (!checkSU2Determinant(m_SU3Generator->generateSU2()))
  {
    cout << "    FAILED: determinant of randomly generated SU2 matrix is not 1."
         << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: random SU2 matrix determinants is 1." << endl;
  }
  return passed;
}

bool SU2Operations::checkSU2Determinant(SU2 H)
{
  /*
   * Function that checks the determinant of the SU2 matrix.
   */
  bool passed = true;
  complex det = SU2Determinant(H);
  if (!((fabs(det.re()) - 1) < m_eps) && !(det.im() < m_eps))
  {
    passed = false;
    if (m_verbose)
      cout << "    FAILED: the determinant of the SU2 matrix differs from 1: "
           << det << endl;
  }
  return passed;
}

bool SU2Operations::run2x2MatrixTests()
{
  /*
   * Function that runs tests on folowing properties of complex 2x2 matrices:
   *  - addition
   *  - subtraction
   *  - multiplication
   *  - conjugation
   *  - matrix transpose
   */
  if (m_verbose && (m_processRank == 0))
    cout << "Running basic SU2 matrix tests." << endl;

  bool passed =
    (testSU2Addition() && testSU2Subtraction() && testSU2Multiplication() &&
     testSU2Conjugation() && testSU2Transpose());

  if (m_processRank == 0)
  {
    if (passed)
    {
      cout << "PASSED: Basic SU2 matrix operations." << endl;
    }
    else
    {
      cout << "FAILURE: Basic SU2 matrix operations." << endl;
    }
  }

  return passed;
}

bool SU2Operations::runSU2Tests()
{
  /*
   * Function that runs tests on folowing properties of SU2 matrices:
   *  - hermicity
   *  - orthogonality of columns
   *  - normality of columns
   *  - determinant equals abs(1)
   */
  bool passed = false;

  if (m_processRank == 0)
  {
    if (m_verbose) cout << "Running SU2 property tests." << endl;

    passed = (testSU2Hermicity() && testSU2Orthogonality() && testSU2Norm() &&
              testSU2Determinant());

    if (passed)
    {
      cout << "PASSED: SU2 properties." << endl;
    }
    else
    {
      cout << "FAILURE: SU2 properties." << endl;
    }
  }

  MPI_Bcast(&passed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
  Parallel::Communicator::setBarrier();

  return passed;
}