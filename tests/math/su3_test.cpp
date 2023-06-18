#include <math/lattice.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <string_view>

const std::string TAGS = "[SU3]";

class SU3Fixture
{
public:
  SU3Fixture()
  {
    // Setting up matrix U1
    mat1.setComplex(complex(1, 1), 0);
    mat1.setComplex(complex(1, 2), 2);
    mat1.setComplex(complex(1, 3), 4);
    mat1.setComplex(complex(2, 1), 6);
    mat1.setComplex(complex(2, 2), 8);
    mat1.setComplex(complex(2, 3), 10);
    mat1.setComplex(complex(3, 1), 12);
    mat1.setComplex(complex(3, 2), 14);
    mat1.setComplex(complex(3, 3), 16);

    // Setting up matrix U2
    mat2.setComplex(complex(4, 4), 0);
    mat2.setComplex(complex(4, 5), 2);
    mat2.setComplex(complex(4, 6), 4);
    mat2.setComplex(complex(5, 4), 6);
    mat2.setComplex(complex(5, 5), 8);
    mat2.setComplex(complex(5, 6), 10);
    mat2.setComplex(complex(6, 4), 12);
    mat2.setComplex(complex(6, 5), 14);
    mat2.setComplex(complex(6, 6), 16);
  }
  ~SU3Fixture() {}

  SU3 mat1, mat2;
};

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
    // Addition
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
  SECTION("Subtraction") {}
  SECTION("Multiplication") {}
  SECTION("Transpose") {}
  SECTION("Conjugation") {}
  SECTION("Complex conjugation") {}

  //   REQUIRE(1 + 1 == 2);
}

TEST_CASE("SU(3) random generator properties")
{
  SECTION("Default random generator")
  {
    SECTION("Hermicity") {}
    SECTION("Orthogonality") {}
    SECTION("Norm") {}
    SECTION("Determinant") {}
  }

  SECTION("RST random generator")
  {
    SECTION("Hermicity") {}
    SECTION("Orthogonality") {}
    SECTION("Norm") {}
    SECTION("Determinant") {}
  }
}

TEST_CASE("SU(3) class methods")
{
  SECTION("Trace") {}
  SECTION("Hermitian") {}
  SECTION("Anti-hermitian") {}
  SECTION("Real trace multiplication") {}
}