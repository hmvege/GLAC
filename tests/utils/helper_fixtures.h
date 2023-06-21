#ifndef GLAC_TESTS_HELPERS_FIXTURES
#define GLAC_TESTS_HELPERS_FIXTURES

#include <math/lattice.h>

class SU3Fixture
{
public:
  SU3Fixture() : mat1(setupMatrix1()), mat2(setupMatrix2()) {}

  ~SU3Fixture() {}

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

  ~SU2Fixture() {}

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
  ~ComplexFixture() {}

  const complex c1, c2;
};

#endif  // GLAC_TESTS_HELPERS_FIXTURES