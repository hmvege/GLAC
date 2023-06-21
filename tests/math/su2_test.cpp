#include <helper_fixtures.h>
#include <helper_generators.h>
#include <helper_math.h>
#include <math/lattice.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <string_view>

const std::string TAGS = "[SU2]";

TEST_CASE_METHOD(SU2Fixture, "SU(2) Matrix tests", TAGS + "[operations]")
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
    SU2 mat_add;
    mat_add.setComplex(complex(5, 5), 0);
    mat_add.setComplex(complex(5, 7), 2);
    mat_add.setComplex(complex(7, 5), 4);
    mat_add.setComplex(complex(7, 7), 6);

    const auto result = mat1 + mat2;
    REQUIRE(result == mat_add);
  }

  SECTION("Subtraction")
  {
    SU2 mat_subtract;
    mat_subtract.setComplex(complex(-3, -3), 0);
    mat_subtract.setComplex(complex(-3, -3), 2);
    mat_subtract.setComplex(complex(-3, -3), 4);
    mat_subtract.setComplex(complex(-3, -3), 6);

    const auto result = mat1 - mat2;
    REQUIRE(result == mat_subtract);
  }

  SECTION("Multiplication")
  {
    SU2 mat_mult;
    mat_mult.setComplex(complex(-3, 22), 0);
    mat_mult.setComplex(complex(-6, 24), 2);
    mat_mult.setComplex(complex(6, 30), 4);
    mat_mult.setComplex(complex(3, 34), 6);

    const auto result = mat1 * mat2;
    REQUIRE(result == mat_mult);
  }

  SECTION("Transpose")
  {
    SU2 mat_transpose;
    mat_transpose.setComplex(complex(1, 1), 0);
    mat_transpose.setComplex(complex(2, 1), 2);
    mat_transpose.setComplex(complex(1, 2), 4);
    mat_transpose.setComplex(complex(2, 2), 6);

    SU2 mat = mat1;

    const auto result = mat.transpose();
    REQUIRE(result == mat_transpose);
  }

  SECTION("Conjugation")
  {
    SU2 mat_conjugate;
    mat_conjugate.setComplex(complex(1, -1), 0);
    mat_conjugate.setComplex(complex(1, -2), 2);
    mat_conjugate.setComplex(complex(2, -1), 4);
    mat_conjugate.setComplex(complex(2, -2), 6);

    SU2 mat = mat1;

    const auto result = mat.conjugate();
    REQUIRE(result == mat_conjugate);
  }

  SECTION("Complex conjugation")
  {
    SU2 mat_complex_conjugate;
    mat_complex_conjugate.setComplex(complex(1, -1), 0);
    mat_complex_conjugate.setComplex(complex(2, -1), 2);
    mat_complex_conjugate.setComplex(complex(1, -2), 4);
    mat_complex_conjugate.setComplex(complex(2, -2), 6);

    SECTION("Complex conjugation then transpose")
    {
      SU2 mat = mat1;
      const auto result = mat.conjugate().transpose();
      REQUIRE(result == mat_complex_conjugate);
    }

    SECTION("Transpose then complex conjugation")
    {
      SU2 mat = mat1;
      const auto result = mat.transpose().conjugate();
      REQUIRE(result == mat_complex_conjugate);
    }
  }
}

TEST_CASE_METHOD(UnitaryMatrixGenerator, "SU(2) random generator properties",
                 TAGS + "[properties]")
{
  const auto mat = generateSU2();

  SECTION("Hermiticity") { REQUIRE(SU2Helpers::checkSU2Hermiticity(mat)); }
  SECTION("Orthogonality") { REQUIRE(SU2Helpers::checkSU2Orthogonality(mat)); }
  SECTION("Norm") { REQUIRE(SU2Helpers::checkSU2Norm(mat)); }
  SECTION("Determinant") { REQUIRE(SU2Helpers::checkSU2Determinant(mat)); }
}