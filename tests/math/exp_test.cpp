#include <helper_fixtures.h>
#include <helper_generators.h>
#include <math/exponentiation/expluscher.h>
#include <math/exponentiation/su3exp.h>
#include <math/exponentiation/taylorexp.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

const std::string TAGS = "[EXPONENTIATION]";

TEST_CASE_METHOD(SU3ExponentiationFixture, "Exponentiation test",
                 TAGS + "[EXP]")
{
  auto exp_method = GENERATE(SU3Exp{}, ExpLuscher{}, TaylorExp{16});

  constexpr double epsilon = 1e-15;

  SECTION("Exponentiation of anti-hermitian matrix")
  {
    // TODO: remove this once we have proper support for constness
    SU3 mat_anti_hermitian = anti_hermitian_matrix;

    const SU3 result_anti_hermitian = exp_method.exp(mat_anti_hermitian);

    for (int i = 0; i < 18; i++)
    {
      REQUIRE_THAT(
        result_anti_hermitian.mat[i],
        Catch::Matchers::WithinAbs(exponentiated_matrix.mat[i], epsilon));
    }
  }

  SECTION(
    "Exponentiation of hermitian matrix and that it yields the same result "
    "as anti-hermitian matrix")
  {
    SU3 mat_hermitian = anti_hermitian_matrix;
    mat_hermitian = mat_hermitian.transpose().conjugate();

    SU3 result_hermitian = exp_method.exp(mat_hermitian);
    SU3 result_transposed_conjugated = result_hermitian;
    result_transposed_conjugated =
      result_transposed_conjugated.transpose().conjugate();

    for (int i = 0; i < 18; i++)
    {
      REQUIRE_THAT(
        result_transposed_conjugated.mat[i],
        Catch::Matchers::WithinAbs(exponentiated_matrix.mat[i], epsilon));
    }
  }
}
