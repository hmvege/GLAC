#include <helper_fixtures.h>
#include <math/matrices/su3.h>
#include <math/matrices/su3matrixgenerator.h>

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <string_view>

const std::string TAGS = "[RSTRandomMatrix]";

class SU3MatrixGeneratorWrapper : public SU3MatrixGenerator
{
public:
  SU3 getRSTMatrixMultiplication(const SU2 &r, const SU2 &s, const SU2 &t)
  {
    return RSTMatrixMultiplication(r, s, t);
  }

  SU3 getRSTMatrixMultiplicationInverse(const SU2 &r, const SU2 &s,
                                        const SU2 &t)
  {
    return RSTMatrixMultiplicationInverse(r, s, t);
  }
};

TEST_CASE_METHOD(MatrixGeneratorFixture, "RST Matrix",
                 TAGS + "[multiplications]")
{
  SU3MatrixGeneratorWrapper generator;

  SECTION("RST matrix multiplication")
  {
    const SU3 results =
      generator.getRSTMatrixMultiplication(mat_r, mat_s, mat_t);

    REQUIRE(results == mat_rst);
  }

  SECTION("RST matrix inverse multiplication")
  {
    const SU3 results =
      generator.getRSTMatrixMultiplicationInverse(mat_r, mat_s, mat_t);

    REQUIRE(results == mat_rst.inv());
  }
}