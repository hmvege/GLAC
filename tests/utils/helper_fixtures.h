#ifndef GLAC_TESTS_HELPERS_FIXTURES
#define GLAC_TESTS_HELPERS_FIXTURES

#include <math/lattice.h>

class SU3Fixture
{
public:
  SU3Fixture() : mat1(setupMatrix1()), mat2(setupMatrix2()) {}

  ~SU3Fixture() = default;

  /** Setting up matrix mat1 */
  SU3 setupMatrix1()
  {
    SU3 temp;
    temp.setComplex(complex(1, 1), 0);
    temp.setComplex(complex(1, 2), 2);
    temp.setComplex(complex(1, 3), 4);
    temp.setComplex(complex(2, 1), 6);
    temp.setComplex(complex(2, 2), 8);
    temp.setComplex(complex(2, 3), 10);
    temp.setComplex(complex(3, 1), 12);
    temp.setComplex(complex(3, 2), 14);
    temp.setComplex(complex(3, 3), 16);

    return temp;
  }

  /** Setting up matrix mat2 */
  SU3 setupMatrix2()
  {
    SU3 temp;
    temp.setComplex(complex(4, 4), 0);
    temp.setComplex(complex(4, 5), 2);
    temp.setComplex(complex(4, 6), 4);
    temp.setComplex(complex(5, 4), 6);
    temp.setComplex(complex(5, 5), 8);
    temp.setComplex(complex(5, 6), 10);
    temp.setComplex(complex(6, 4), 12);
    temp.setComplex(complex(6, 5), 14);
    temp.setComplex(complex(6, 6), 16);

    return temp;
  }

  const SU3 mat1, mat2;
};

class SU2Fixture
{
public:
  SU2Fixture() : mat1(setupMatrix1()), mat2(setupMatrix2()) {}

  ~SU2Fixture() = default;

  /** Setting up matrix mat1 */
  SU2 setupMatrix1()
  {
    SU2 temp;
    temp.setComplex(complex(1, 1), 0);
    temp.setComplex(complex(1, 2), 2);
    temp.setComplex(complex(2, 1), 4);
    temp.setComplex(complex(2, 2), 6);

    return temp;
  }

  /** Setting up matrix mat2 */
  SU2 setupMatrix2()
  {
    SU2 temp;
    temp.setComplex(complex(4, 4), 0);
    temp.setComplex(complex(4, 5), 2);
    temp.setComplex(complex(5, 4), 4);
    temp.setComplex(complex(5, 5), 6);

    return temp;
  }

  const SU2 mat1, mat2;
};

class ComplexFixture
{
public:
  ComplexFixture() : c1(complex(1, 2)), c2(complex(3, 4)) {}
  ~ComplexFixture() = default;

  const complex c1, c2;
};

class MatrixGeneratorFixture
{
public:
  MatrixGeneratorFixture() :
    mat_rst(setupRST()), mat_r(setupR()), mat_s(setupS()), mat_t(setupT())
  {}
  ~MatrixGeneratorFixture() = default;

  const SU3 mat_rst;
  const SU2 mat_r, mat_s, mat_t;

private:
  SU3 setupRST()
  {
    SU3 matrix_rst;
    matrix_rst.setComplex({-8, 9}, 0);
    matrix_rst.setComplex({-123, 82}, 2);
    matrix_rst.setComplex({-78, 83}, 4);
    matrix_rst.setComplex({-20, 37}, 6);
    matrix_rst.setComplex({-391, 342}, 8);
    matrix_rst.setComplex({-234, 319}, 10);
    matrix_rst.setComplex({2, 6}, 12);
    matrix_rst.setComplex({48, 58}, 14);
    matrix_rst.setComplex({44, 35}, 16);
    return matrix_rst;
  }

  SU2 setupR()
  {
    SU2 mat;
    mat.setComplex(complex(1, 2), 0);
    mat.setComplex(complex(3, 4), 2);
    mat.setComplex(complex(5, 6), 4);
    mat.setComplex(complex(7, 8), 6);
    return mat;
  }

  SU2 setupS()
  {
    SU2 mat;
    mat.setComplex(complex(2, 5), 0);
    mat.setComplex(complex(4, 8), 2);
    mat.setComplex(complex(2, 6), 4);
    mat.setComplex(complex(10, 3), 6);
    return mat;
  }

  SU2 setupT()
  {
    SU2 mat;
    mat.setComplex(complex(7, 2), 0);
    mat.setComplex(complex(6, 1), 2);
    mat.setComplex(complex(6, 4), 4);
    mat.setComplex(complex(5, 2), 6);
    return mat;
  }
};

#endif  // GLAC_TESTS_HELPERS_FIXTURES