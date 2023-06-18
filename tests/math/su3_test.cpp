#include <helper_fixtures.h>
#include <helper_generators.h>
#include <helper_math.h>
#include <math/lattice.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <string_view>

const std::string TAGS = "[SU3]";

TEST_CASE_METHOD(SU3Fixture, "SU(3) Matrix tests", TAGS + "[operations]")
{
  SECTION("Equality")
  {
    REQUIRE(mat1 == mat1);
    REQUIRE_FALSE(mat1 != mat1);

    REQUIRE(mat2 == mat2);
    REQUIRE_FALSE(mat2 != mat2);

    REQUIRE(mat1 != mat2);
    REQUIRE_FALSE(mat1 == mat2);
  }

  SECTION("Addition")
  {
    SU3 mat_add;
    mat_add.setComplex(complex(5, 5), 0);
    mat_add.setComplex(complex(5, 7), 2);
    mat_add.setComplex(complex(5, 9), 4);
    mat_add.setComplex(complex(7, 5), 6);
    mat_add.setComplex(complex(7, 7), 8);
    mat_add.setComplex(complex(7, 9), 10);
    mat_add.setComplex(complex(9, 5), 12);
    mat_add.setComplex(complex(9, 7), 14);
    mat_add.setComplex(complex(9, 9), 16);

    const auto result = mat1 + mat2;
    REQUIRE(result == mat_add);
  }

  SECTION("Subtraction")
  {
    SU3 mat_subtract;
    mat_subtract.setComplex(complex(-3, -3), 0);
    mat_subtract.setComplex(complex(-3, -3), 6);
    mat_subtract.setComplex(complex(-3, -3), 12);
    mat_subtract.setComplex(complex(-3, -3), 2);
    mat_subtract.setComplex(complex(-3, -3), 8);
    mat_subtract.setComplex(complex(-3, -3), 14);
    mat_subtract.setComplex(complex(-3, -3), 4);
    mat_subtract.setComplex(complex(-3, -3), 10);
    mat_subtract.setComplex(complex(-3, -3), 16);

    const auto result = mat1 - mat2;
    REQUIRE(result == mat_subtract);
  }

  SECTION("Multiplication")
  {
    SU3 mat_mult;
    mat_mult.setComplex(complex(-9, 44), 0);
    mat_mult.setComplex(complex(-15, 47), 2);
    mat_mult.setComplex(complex(-21, 50), 4);
    mat_mult.setComplex(complex(6, 56), 6);
    mat_mult.setComplex(complex(0, 62), 8);
    mat_mult.setComplex(complex(-6, 68), 10);
    mat_mult.setComplex(complex(21, 68), 12);
    mat_mult.setComplex(complex(15, 77), 14);
    mat_mult.setComplex(complex(9, 86), 16);

    const auto result = mat1 * mat2;
    REQUIRE(result == mat_mult);
  }

  SECTION("Transpose")
  {
    SU3 mat_transpose;
    mat_transpose.setComplex(complex(1, 1), 0);
    mat_transpose.setComplex(complex(2, 1), 2);
    mat_transpose.setComplex(complex(3, 1), 4);
    mat_transpose.setComplex(complex(1, 2), 6);
    mat_transpose.setComplex(complex(2, 2), 8);
    mat_transpose.setComplex(complex(3, 2), 10);
    mat_transpose.setComplex(complex(1, 3), 12);
    mat_transpose.setComplex(complex(2, 3), 14);
    mat_transpose.setComplex(complex(3, 3), 16);

    SU3 mat = mat1;

    const auto result = mat.transpose();
    REQUIRE(result == mat_transpose);
  }

  SECTION("Conjugation")
  {
    SU3 mat_conjugate;
    mat_conjugate.setComplex(complex(1, -1), 0);
    mat_conjugate.setComplex(complex(1, -2), 2);
    mat_conjugate.setComplex(complex(1, -3), 4);
    mat_conjugate.setComplex(complex(2, -1), 6);
    mat_conjugate.setComplex(complex(2, -2), 8);
    mat_conjugate.setComplex(complex(2, -3), 10);
    mat_conjugate.setComplex(complex(3, -1), 12);
    mat_conjugate.setComplex(complex(3, -2), 14);
    mat_conjugate.setComplex(complex(3, -3), 16);

    SU3 mat = mat1;

    const auto result = mat.conjugate();
    REQUIRE(result == mat_conjugate);
  }

  SECTION("Complex conjugation")
  {
    SU3 mat_complex_conjugate;
    mat_complex_conjugate.setComplex(complex(1, -1), 0);
    mat_complex_conjugate.setComplex(complex(2, -1), 2);
    mat_complex_conjugate.setComplex(complex(3, -1), 4);
    mat_complex_conjugate.setComplex(complex(1, -2), 6);
    mat_complex_conjugate.setComplex(complex(2, -2), 8);
    mat_complex_conjugate.setComplex(complex(3, -2), 10);
    mat_complex_conjugate.setComplex(complex(1, -3), 12);
    mat_complex_conjugate.setComplex(complex(2, -3), 14);
    mat_complex_conjugate.setComplex(complex(3, -3), 16);

    SECTION("Complex conjugation then transpose")
    {
      SU3 mat = mat1;
      const auto result = mat.conjugate().transpose();
      REQUIRE(result == mat_complex_conjugate);
    }

    SECTION("Transpose then complex conjugation")
    {
      SU3 mat = mat1;
      const auto result = mat.transpose().conjugate();
      REQUIRE(result == mat_complex_conjugate);
    }
  }
}

TEST_CASE_METHOD(UnitaryMatrixGenerator, "SU(3) random generator properties",
                 TAGS + "[properties]")
{
  SECTION("Default random generator")
  {
    const auto mat = generator.generateRandom();

    SECTION("Hermiticity") { REQUIRE(checkSU3Hermiticity(mat)); }
    SECTION("Orthogonality") { REQUIRE(checkSU3Orthogonality(mat)); }
    SECTION("Norm") { REQUIRE(checkSU3Norm(mat)); }
    SECTION("Determinant") { REQUIRE(checkSU3Determinant(mat)); }
  }

  SECTION("RST random generator")
  {
    const auto mat = generator.generateRST();

    SECTION("Hermiticity") { REQUIRE(checkSU3Hermiticity(mat)); }
    SECTION("Orthogonality") { REQUIRE(checkSU3Orthogonality(mat)); }
    SECTION("Norm") { REQUIRE(checkSU3Norm(mat)); }
    SECTION("Determinant") { REQUIRE(checkSU3Determinant(mat)); }
  }
}

TEST_CASE_METHOD(SU3Fixture, "SU(3) class methods", TAGS + "[su3methods]")
{
  SECTION("Trace")
  {
    complex mat1_trace(6, 6);

    SU3 mat = mat1;

    REQUIRE(mat.trace() == mat1_trace);
  }

  SECTION("Hermiticity")
  {
    SU3 mat_hermitian;
    mat_hermitian.setComplex(complex(-4, 0), 0);
    mat_hermitian.setComplex(complex(-3, 1), 2);
    mat_hermitian.setComplex(complex(-2, 2), 4);
    mat_hermitian.setComplex(complex(-3, -1), 6);
    mat_hermitian.setComplex(complex(-2, 0), 8);
    mat_hermitian.setComplex(complex(-1, 1), 10);
    mat_hermitian.setComplex(complex(-2, -2), 12);
    mat_hermitian.setComplex(complex(-1, -1), 14);
    mat_hermitian.setComplex(complex(0, 0), 16);

    SU3 mat_anti_hermitian;
    mat_anti_hermitian.setComplex(complex(0, -4), 0);
    mat_anti_hermitian.setComplex(complex(-1, -3), 2);
    mat_anti_hermitian.setComplex(complex(-2, -2), 4);
    mat_anti_hermitian.setComplex(complex(1, -3), 6);
    mat_anti_hermitian.setComplex(complex(0, -2), 8);
    mat_anti_hermitian.setComplex(complex(-1, -1), 10);
    mat_anti_hermitian.setComplex(complex(2, -2), 12);
    mat_anti_hermitian.setComplex(complex(1, -1), 14);
    mat_anti_hermitian.setComplex(complex(0, 0), 16);

    SECTION("Hermitian")
    {
      SU3 mat = mat_anti_hermitian;
      mat.makeHermitian();

      REQUIRE(is_close_su3(mat, mat_hermitian));
    }

    SECTION("Anti-hermitian")
    {
      SU3 mat = mat_hermitian;
      mat.makeAntiHermitian();

      REQUIRE(is_close_su3(mat, mat_anti_hermitian));
    }
  }

  SECTION("Real trace multiplication")
  {
    double traced_matrix = -21.0;  // re(tr(U1*U3))

    SU3 mat_trace;
    mat_trace.setComplex(complex(4, 4), 0);
    mat_trace.setComplex(complex(5, 5), 2);
    mat_trace.setComplex(complex(2, 6), 4);
    mat_trace.setComplex(complex(6, 4), 6);
    mat_trace.setComplex(complex(5, 5), 8);
    mat_trace.setComplex(complex(1, 6), 10);
    mat_trace.setComplex(complex(6, 4), 12);
    mat_trace.setComplex(complex(9, 5), 14);
    mat_trace.setComplex(complex(2, 6), 16);

    SU3 mat = mat1;

    const auto results = traceRealMultiplication(mat, mat_trace);

    REQUIRE(traced_matrix == results);
  }
}