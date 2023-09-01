#include <helper_fixtures.h>
#include <helper_generators.h>
#include <helper_math.h>
#include <math/functions.h>

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <random>
#include <stdexcept>

const std::string TAGS = "[FUNCTIONS]";

TEST_CASE_METHOD(DeterminantMatrices, "Determinant verification",
                 TAGS + "[DETERMINANTS]")
{
  constexpr double epsilon = 1e-14;

  SECTION("2x2 determinant")
  {
    SECTION("SU2 matrix")
    {
      const complex determinant = SU2Determinant(su2_matrix);
      CAPTURE(determinant);
      REQUIRE(
        ComplexHelpers::isWithinAbs(determinant, complex(1.0, 0.0), epsilon));
    }
    SECTION("Random matrix")
    {
      const complex determinant = SU2Determinant(random_2x2_matrix);
      CAPTURE(determinant);
      REQUIRE(ComplexHelpers::isWithinAbs(
        determinant, complex(0.10144648144193713, -1.6677931532846635),
        epsilon));
    }
  }

  SECTION("3x3 determinant")
  {
    SECTION("SU3 matrix")
    {
      const complex determinant = SU3Determinant(su3_matrix);
      CAPTURE(determinant);
      REQUIRE(
        ComplexHelpers::isWithinAbs(determinant, complex(1.0, 0.0), epsilon));
    }
    SECTION("Random matrix")
    {
      const complex determinant = SU3Determinant(random_3x3_matrix);
      CAPTURE(determinant);
      REQUIRE(ComplexHelpers::isWithinAbs(
        determinant, complex(-8.406173922755011, 1.238596912388841), epsilon));
    }
  }
}

TEST_CASE_METHOD(DeterminantMatrices, "Trace multiplication verification",
                 TAGS + "[TRACE]")
{
  constexpr double epsilon = 1e-15;

  SECTION("Real trace of two SU3 matrices, A * B")
  {
    constexpr double real_trace = -1.053769558618308;

    SECTION("Testing with A and B as reference")
    {
      const SU3 matrix_a = su3_matrix;
      const SU3 matrix_b = random_3x3_matrix;
      REQUIRE_THAT(traceRealMultiplication(matrix_a, matrix_b),
                   Catch::Matchers::WithinAbs(real_trace, epsilon));
    }
    SECTION("Testing with A as an r-value and B as reference")
    {
      SU3 matrix_a = su3_matrix;
      const SU3 matrix_b = random_3x3_matrix;
      REQUIRE_THAT(traceRealMultiplication(std::move(matrix_a), matrix_b),
                   Catch::Matchers::WithinAbs(real_trace, epsilon));
    }
    SECTION("Testing with A as reference and B as an r-value reference")
    {
      const SU3 matrix_a = su3_matrix;
      SU3 matrix_b = random_3x3_matrix;
      REQUIRE_THAT(traceRealMultiplication(matrix_a, std::move(matrix_b)),
                   Catch::Matchers::WithinAbs(real_trace, epsilon));
    }
  }

  SECTION("Imaginary trace of two SU3 matrices, A * B")
  {
    constexpr double imag_trace = 2.501934167586923;

    SECTION("Testing with A and B as reference")
    {
      const SU3 matrix_a = su3_matrix;
      const SU3 matrix_b = random_3x3_matrix;
      REQUIRE_THAT(traceImagMultiplication(matrix_a, matrix_b),
                   Catch::Matchers::WithinAbs(imag_trace, epsilon));
    }
    SECTION("Testing with A as an r-value and B as reference")
    {
      SU3 matrix_a = su3_matrix;
      const SU3 matrix_b = random_3x3_matrix;
      REQUIRE_THAT(traceImagMultiplication(std::move(matrix_a), matrix_b),
                   Catch::Matchers::WithinAbs(imag_trace, epsilon));
    }
    SECTION("Testing with A as reference and B as an r-value reference")
    {
      const SU3 matrix_a = su3_matrix;
      SU3 matrix_b = random_3x3_matrix;
      REQUIRE_THAT(traceImagMultiplication(matrix_a, std::move(matrix_b)),
                   Catch::Matchers::WithinAbs(imag_trace, epsilon));
    }
  }
}
