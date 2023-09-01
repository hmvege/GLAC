#include <helper_math.h>
#include <math/functions.h>
#include <math/matrices/su3.h>

#include <catch2/catch_all.hpp>
#include <cmath>

namespace ComplexHelpers
{
  bool isWithinAbs(const complex& a, const complex& b, double epsilon)
  {
    return (std::abs(a.re() - b.re()) <= epsilon) &&
           (std::abs(a.im() - b.im()) <= epsilon);
  }

  // Inline dot product between two columns of a 3x3 matrix
  inline complex dot(complex* a, complex* b)
  {
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex return_sum(0, 0);
    for (int i = 0; i < 3; i++)
    {
      return_sum += a[i].c() * b[i];
    }
    return return_sum;
  }

  // Inline dot product between two columns of a 2x2 matrix
  inline complex dot2(complex* a, complex* b)
  {
    /*
     * Dot product defined as:
     *  (v,u*) = v*conjugate(u)
     */
    complex return_sum(0, 0);
    for (int i = 0; i < 2; i++)
    {
      return_sum += a[i].c() * b[i];
    }
    return return_sum;
  }
}  // namespace ComplexHelpers

namespace SU3Helpers
{
  bool checkSU3Hermiticity(SU3 H)
  {
    /*
     * Checks the hermiticity of the SU3 matrices.
     */
    constexpr double eps = 2 * 1e-14;

    bool passed = true;

    SU3 I;
    I = H * H.inv();
    for (int i = 0; i < 18; i++)
    {
      if (i == 0 || i == 8 || i == 16)
      {
        continue;
      }
      if (std::fabs(I[i]) > eps)
      {
        passed = false;
      }
    }

    if ((fabs(I[0] - 1) > eps) || (fabs(I[8] - 1) > eps) ||
        (fabs(I[16] - 1) > eps))
    {
      passed = false;
    }

    return passed;
  }

  bool checkSU3Orthogonality(SU3 H)
  {
    /*
     * Testing the orthogonality of a SU3 matrix H.
     */
    constexpr double eps = 2 * 1e-14;

    bool passed = true;

    complex col1[3];
    complex col2[3];
    complex col3[3];
    for (int i = 0; i < 3; i++)
    {
      col1[i].setRe(H[6 * i]);
      col1[i].setIm(H[6 * i + 1]);
      col2[i].setRe(H[6 * i + 2]);
      col2[i].setIm(H[6 * i + 3]);
      col3[i].setRe(H[6 * i + 4]);
      col3[i].setIm(H[6 * i + 5]);
    }
    complex c12dot = ComplexHelpers::dot(col1, col2);
    complex c13dot = ComplexHelpers::dot(col1, col3);
    complex c23dot = ComplexHelpers::dot(col2, col3);
    if ((std::fabs(c12dot.re()) > eps) || (fabs(c12dot.im()) > eps))
    {
      passed = false;
    }
    if ((std::fabs(c13dot.re()) > eps) || (fabs(c13dot.im()) > eps))
    {
      passed = false;
    }
    if ((std::fabs(c23dot.re()) > eps) || (fabs(c23dot.im()) > eps))
    {
      passed = false;
    }

    return passed;
  }

  bool checkSU3Norm(SU3 H)
  {
    /*
     * Function that checks the SU3 norm of a matrix H.
     */
    constexpr double eps = 2 * 1e-14;

    bool passed = true;

    double sum;
    for (int col = 0; col < 3; col++)
    {
      sum = 0;
      for (int i = 0; i < 3; i++)
      {
        sum += H.normSquared(6 * i + 2 * col);
      }
      sum = sqrt(sum);
      if (std::fabs(sqrt(sum) - 1) > eps)
      {
        passed = false;
      }
    }
    return passed;
  }

  bool checkSU3Determinant(SU3 H)
  {
    /*
     * Function that checks the determinant of the SU2 matrix.
     */
    constexpr double eps = 2 * 1e-14;

    bool passed = true;

    complex det = SU3Determinant(H);

    if (!((std::fabs(det.re()) - 1) < eps) && !(det.im() < eps))
    {
      passed = false;
    }

    return passed;
  }

  bool is_close_su3(const SU3& A, const SU3& B, const double eps)
  {
    bool is_close = true;
    for (int i = 0; i < 18; i++)
    {
      if (std::fabs(A.mat[i] - B.mat[i]) > eps)
      {
        is_close = false;
        break;
      }
    }
    return is_close;
  }
}  // namespace SU3Helpers

namespace SU2Helpers
{
  bool checkSU2Hermiticity(SU2 H)
  {
    /*
     * Checks the hermiticity of the SU2 matrices.
     */
    constexpr double eps = 2 * 1e-14;

    bool passed = true;

    SU2 I;
    I = H * H.inv();
    for (int i = 0; i < 8; i++)
    {
      if (i == 0 || i == 6)
      {
        continue;
      }
      if (std::fabs(I[i]) > eps)
      {
        passed = false;
      }
    }
    if ((std::fabs(I[0] - 1) > eps) || (std::fabs(I[6] - 1) > eps))
    {
      passed = false;
    }

    return passed;
  }

  bool checkSU2Orthogonality(SU2 H)
  {
    /*
     * Testing the orthogonality of a SU3 matrix H.
     */
    constexpr double eps = 2 * 1e-14;

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
    complex c12dot = ComplexHelpers::dot2(col1, col2);
    if ((std::fabs(c12dot.re()) > eps) || (std::fabs(c12dot.im()) > eps))
    {
      passed = false;
    }

    return passed;
  }

  bool checkSU2Norm(SU2 H)
  {
    /*
     * Function that checks the SU2 norm of a matrix H.
     */
    constexpr double eps = 2 * 1e-14;

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
      if (std::fabs(sqrt(sum) - 1) > eps)
      {
        passed = false;
      }
    }

    return passed;
  }

  bool checkSU2Determinant(SU2 H)
  {
    /*
     * Function that checks the determinant of the SU2 matrix.
     */
    constexpr double eps = 2 * 1e-14;

    bool passed = true;

    complex det = SU2Determinant(H);

    if (!((std::fabs(det.re()) - 1) < eps) && !(det.im() < eps))
    {
      passed = false;
    }

    return passed;
  }

  bool is_close_su2(const SU2& A, const SU2& B, const double eps)
  {
    bool is_close = true;
    for (int i = 0; i < 8; i++)
    {
      if (std::fabs(A.mat[i] - B.mat[i]) > eps)
      {
        is_close = false;
        break;
      }
    }
    return is_close;
  }
}  // namespace SU2Helpers