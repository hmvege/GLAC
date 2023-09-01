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
    // TODO: implement a direct constructor for SU3
    // Currently, this method is rather unfortunate, since it would be better
    // to implement a constructor which can directly initialise matrices from
    // e.g. vector of vectors of complex numbers or similar.
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

/**
 * \brief Helper fixture for the SU3 exponentiation function
 *
 * The SU3 exponentiation requires a traceless, hermitian matrix.
 */
class SU3ExponentiationFixture
{
public:
  SU3ExponentiationFixture() :
    anti_hermitian_matrix(setupAntiHermitianMatrix()),
    exponentiated_matrix(setupExponentiatedMatrix())
  {}
  ~SU3ExponentiationFixture() = default;

  /** An SU3 RST-generated matrix */
  const SU3 anti_hermitian_matrix;
  /** The resulting exponentiated \p anti_hermitian_matrix */
  const SU3 exponentiated_matrix;

private:
  SU3 setupAntiHermitianMatrix()
  {
    SU3 mat;
    mat.setComplex({0.0000000000000000, -0.0973617343761928}, 0);
    mat.setComplex({0.1940344580490270, -0.1826110044731708}, 2);
    mat.setComplex({0.0006214121141749, -0.0014679912989845}, 4);
    mat.setComplex({-0.1940344580490270, -0.1826110044731708}, 6);
    mat.setComplex({0.0000000000000000, 0.0793466046243053}, 8);
    mat.setComplex({-0.0314735262421320, -0.0131822500321111}, 10);
    mat.setComplex({-0.0006214121141749, -0.0014679912989845}, 12);
    mat.setComplex({0.0314735262421320, -0.0131822500321111}, 14);
    mat.setComplex({0.0000000000000000, 0.0180151297518875}, 16);
    return mat;
  }

  SU3 setupExponentiatedMatrix()
  {
    SU3 mat;
    mat.setComplex({0.9600346049322120, -0.0958460957487774}, 0);
    mat.setComplex({0.1898270929729914, -0.1819319982836562}, 2);
    mat.setComplex({-0.0036733740024759, 0.0001111743883831}, 4);
    mat.setComplex({-0.1930746457073775, -0.1784270410252529}, 6);
    mat.setComplex({0.9610289600574021, 0.0785093402600774}, 8);
    mat.setComplex({-0.0306092729886854, -0.0144452776353779}, 10);
    mat.setComplex({-0.0048995805100383, -0.0030082535027400}, 12);
    mat.setComplex({0.0314995566655486, -0.0115681513963342}, 14);
    mat.setComplex({0.9992582585838492, 0.0179953598687933}, 16);
    return mat;
  }
};

/**
 * \brief Helper fixture for calculating and verifying determinants
 */
class DeterminantMatrices
{
public:
  DeterminantMatrices() :
    su3_matrix(setupSu3Matrix()),
    su2_matrix(setupSu2Matrix()),
    random_3x3_matrix(setupRandom3x3Matrix()),
    random_2x2_matrix(setupRandom2x2Matrix())
  {}
  ~DeterminantMatrices() = default;

  /** Sample SU3 matrix for testing the verification */
  const SU3 su3_matrix;
  /** Sample SU2 matrix for testing the verification */
  const SU2 su2_matrix;
  /** Sample random 3x3 matrix for testing the verification */
  const SU3 random_3x3_matrix;
  /** Sample random 2x2 matrix for testing the verification */
  const SU2 random_2x2_matrix;

private:
  SU3 setupSu3Matrix()
  {
    SU3 mat;
    mat.setComplex({9.4570982496843226e-01, 4.3990300027111753e-02}, 0);
    mat.setComplex({7.6168903553637285e-02, -1.7636031344099676e-01}, 2);
    mat.setComplex({-8.2725980109240357e-02, 2.4484593247815559e-01}, 4);
    mat.setComplex({-4.6268049132209228e-02, -2.2604309093565370e-01}, 6);
    mat.setComplex({9.5401180274076880e-01, -1.3308930482705782e-02}, 8);
    mat.setComplex({-1.3672682376043147e-01, -1.3324382512842381e-01}, 10);
    mat.setComplex({4.1749668907258682e-02, 2.2072275432676039e-01}, 12);
    mat.setComplex({1.8261576581539712e-01, -1.3939574115015468e-01}, 14);
    mat.setComplex({9.4619301577331516e-01, -3.8438485859839472e-02}, 16);
    return mat;
  }

  SU3 setupRandom3x3Matrix()
  {
    SU3 mat;
    mat.setComplex({1.8832289713054970e-01, 1.3715556962794211e+00}, 0);
    mat.setComplex({-7.0230687282911686e-01, 4.8366810289453366e-01}, 2);
    mat.setComplex({9.4015870445629057e-01, -2.1619589188205972e+00}, 4);
    mat.setComplex({1.5114463350850751e+00, 1.5856946998978356e+00}, 6);
    mat.setComplex({-8.7250317417400425e-01, -9.7197178917090427e-02}, 8);
    mat.setComplex({3.1276045466995661e-01, 6.6560436959509706e-01}, 10);
    mat.setComplex({-9.4478001524210620e-01, -7.9082448979219977e-01}, 12);
    mat.setComplex({-3.2832084478126244e-01, -8.1801631123549468e-01}, 14);
    mat.setComplex({-1.8921193407330390e+00, 1.0879946448651587e+00}, 16);
    return mat;
  }

  SU2 setupSu2Matrix()
  {
    SU2 mat;
    mat.setComplex({-9.7077288796092776e-01, 3.4361286417098749e-02}, 0);
    mat.setComplex({-1.6227555581986752e-01, 1.7345300798520380e-01}, 2);
    mat.setComplex({1.6227555581986752e-01, 1.7345300798520380e-01}, 4);
    mat.setComplex({-9.7077288796092776e-01, -3.4361286417098749e-02}, 6);
    return mat;
  }

  SU2 setupRandom2x2Matrix()
  {
    SU2 mat;
    mat.setComplex({1.4478863174834373e+00, -1.4831561936819222e+00}, 0);
    mat.setComplex({-6.9610406417799953e-01, 6.5077391558030406e-01}, 2);
    mat.setComplex({1.4720312169753980e-01, 3.9435368177029045e-01}, 4);
    mat.setComplex({5.5063592133706540e-01, -7.1126393858800196e-01}, 6);
    return mat;
  }
};

#endif  // GLAC_TESTS_HELPERS_FIXTURES