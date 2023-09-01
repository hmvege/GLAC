#include <helper_fixtures.h>
#include <helper_generators.h>
#include <math/exponentiation/taylorexp.h>

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <stdexcept>

const std::string TAGS = "[EXPONENTIATION]";

TEST_CASE_METHOD(SU3ExponentiationFixture, "Exponentiation test",
                 TAGS + "[TaylorExpMethods]")
{
  SECTION("Degree limitations")
  {
    // Test if the degree limitations are correctly enforced.
    REQUIRE_THROWS_MATCHES(
      TaylorExp(0), std::invalid_argument,
      Catch::Matchers::Message("Error: invalid degree 0 of Taylor expansion. "
                               "Minimum degree allowed is 1."));
    REQUIRE_THROWS_MATCHES(
      TaylorExp(17), std::invalid_argument,
      Catch::Matchers::Message("Error: invalid degree 17 of Taylor expansion. "
                               "Maximum allowed degree is 16"));
  }
}
