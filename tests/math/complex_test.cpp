#include <helper_fixtures.h>
// #include <helper_generators.h>
// #include <helper_math.h>
#include <math/complex.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <string_view>

const std::string TAGS = "[Complex]";

TEST_CASE_METHOD(ComplexFixture, "Complex unit tests", TAGS + "[operations]")
{
  SECTION("Equality")
  {
    REQUIRE(c1 == c1);
    REQUIRE_FALSE(c1 != c1);

    REQUIRE(c2 == c2);
    REQUIRE_FALSE(c2 != c2);

    REQUIRE(c1 != c2);
    REQUIRE_FALSE(c1 == c2);
  }

  SECTION("Addition")
  {
    const complex c_add(4, 6);
    REQUIRE((c1 + c2) == c_add);
  }

  SECTION("Subtraction")
  {
    const complex c_sub(-2, -2);
    REQUIRE((c1 - c2) == c_sub);
  }

  SECTION("Multiplication")
  {
    const complex c_mul(-5, 10);
    REQUIRE(c1 * c2 == c_mul);
  }

  SECTION("Division")
  {
    const complex c_div(0.44, 0.08);
    REQUIRE(c1 / c2 == c_div);
  }

  SECTION("Set to minus operator")
  {
    const complex c_set_to_minus(-1, -2);
    complex temp;
    temp = -c1;
    REQUIRE(temp == c_set_to_minus);
  }

  SECTION("Multiplying by 0-i")
  {
    const complex c_minus_i(0, -1);
    const complex expected(2, -1);
    REQUIRE(c1 * c_minus_i == expected);
  }
}

TEST_CASE_METHOD(ComplexFixture, "Complex properties", TAGS + "[properties]")
{
  SECTION("Conjugate")
  {
    const complex c_conj(1, -2);
    REQUIRE(c1.c() == c_conj);

    auto temp = c1;
    REQUIRE(temp.conjugate() == c_conj);
  }

  SECTION("Norm")
  {
    const double c_norm = 2.23606797749979;
    REQUIRE(c1.norm() == c_norm);
  }

  SECTION("Norm Squared")
  {
    const double c_norm_squared = 5;
    REQUIRE(c1.normSquared() == c_norm_squared);
  }
}