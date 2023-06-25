#include "complexoperations.h"

#include "parallelization/communicator.h"

ComplexOperations::ComplexOperations() {}

////////////////////////////////////
////////// COMPLEX TESTS ///////////
////////////////////////////////////
bool ComplexOperations::testComplexAddition()
{
  bool passed = true;
  if (!compareComplex(z1 + z2, zAdd))
  {
    cout << "    FAILED: complex addition is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex addition is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexSubtraction()
{
  bool passed = true;
  if (!compareComplex(z1 - z2, zSub))
  {
    cout << "    FAILED: complex subtraction is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex subtraction is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexMultiplication()
{
  bool passed = true;
  if (!compareComplex(z1 * z2, zMul))
  {
    cout << "    FAILED: complex multipliction is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex multiplication is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexDivision()
{
  bool passed = true;
  if (!compareComplex(z1 / z2, zDiv))
  {
    cout << "    FAILED: complex division is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex division is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexConjugation()
{
  bool passed = true;
  complex temp = z1;
  if (!compareComplex(temp.c(), zConj))
  {
    cout << "    FAILED: complex conjugation c() is not correct." << endl;
    passed = false;
    temp = z1;
  }
  else if (!compareComplex(temp.conjugate(), zConj))
  {
    cout << "    FAILED: complex conjugation conjugate() is not correct."
         << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex conjugate() and c() is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexNorm()
{
  bool passed = true;
  if (fabs(z1.norm() - zNorm) > m_machineEpsilon)
  {
    cout << "    FAILED: complex norm is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex norm is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexNormSquared()
{
  bool passed = true;
  if (fabs(z1.normSquared() - zNormSquared) > m_machineEpsilon)
  {
    cout << "    FAILED: complex normSquared is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: complex normSquared is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::testComplexSetToMinus()
{
  bool passed = true;
  complex temp;
  temp = -z1;
  if (!compareComplex(temp, zSetToMinus))
  {
    cout << "    FAILED: setting complex to minus is not correct." << endl;
    passed = false;
  }
  if (passed && m_verbose)
  {
    cout << "    SUCCESS: setting complex to minus is correct." << endl;
  }
  return passed;
}

bool ComplexOperations::runComplexTests()
{
  /*
   * Function for testing all aspects of the complex class:
   *  - addition
   *  - subtraction
   *  - multiplication
   *  - division
   *  - conjugation with conjugate()
   *  - conjugation with c()
   *  - set to minus operator
   *  - norm
   *  - norm squared
   */
  if (m_verbose && (m_processRank == 0))
    cout << "Running complex tests." << endl;
  bool passed = false;

  if (m_processRank == 0)
  {
    passed = (testComplexAddition() && testComplexSubtraction() &&
              testComplexMultiplication() && testComplexDivision() &&
              testComplexConjugation() && testComplexNorm() &&
              testComplexNormSquared() && testComplexSetToMinus());
    if (passed)
    {
      cout << "PASSED: Complex operations." << endl;
    }
    else
    {
      cout << "FAILURE: Complex operations." << endl;
    }
  }
  MPI_Bcast(&passed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
  Parallel::Communicator::setBarrier();

  return passed;
}
