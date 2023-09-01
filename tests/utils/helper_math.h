#ifndef GLAC_TESTS_HELPERS_MATH
#define GLAC_TESTS_HELPERS_MATH

#include <math/matrices/su2.h>
#include <math/matrices/su3.h>

namespace SU3Helpers
{
  bool checkSU3Hermiticity(SU3 H);

  bool checkSU3Orthogonality(SU3 H);

  bool checkSU3Norm(SU3 H);

  bool checkSU3Determinant(SU3 H);

  bool is_close_su3(const SU3& A, const SU3& B, const double eps = 1e-15);
}  // namespace SU3Helpers

namespace SU2Helpers
{
  bool checkSU2Hermiticity(SU2 H);

  bool checkSU2Orthogonality(SU2 H);

  bool checkSU2Norm(SU2 H);

  bool checkSU2Determinant(SU2 H);

  bool is_close_su2(const SU2& A, const SU2& B, const double eps = 1e-15);

}  // namespace SU2Helpers

namespace ComplexHelpers
{
  /** Checks that two complex numbers are within an epsilon value of each
   * other*/
  bool isWithinAbs(const complex& a, const complex& b, double epsilon);

  /** Inline dot product between two columns of a 3x3 matrix */
  inline complex dot(complex* a, complex* b);
  /** Inline dot product between two columns of a 2x2 matrix */
  inline complex dot2(complex* a, complex* b);

}  // namespace ComplexHelpers
#endif  // GLAC_TESTS_HELPERS_MATH