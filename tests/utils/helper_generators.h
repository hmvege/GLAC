#ifndef GLAC_TESTS_HELPERS_GENERATOR
#define GLAC_TESTS_HELPERS_GENERATOR

#include <math/matrices/su3matrixgenerator.h>

class UnitaryMatrixGenerator
{
public:
  UnitaryMatrixGenerator() { generator = SU3MatrixGenerator(); }
  ~UnitaryMatrixGenerator() = default;

  SU2 generateSU2() { return generator.generateSU2(); }
  SU3 generateSU3() { return generator.generateRandom(); }
  SU3 generateSU3RST() { return generator.generateRST(); }

private:
  SU3MatrixGenerator generator;
};

class RandomGenerator
{
public:
  RandomGenerator()
  {
    m_generator = std::mt19937_64(0);
    m_uniform_distribution = std::uniform_real_distribution<double>(
      0, 1);  // Print out max values to ensure we dont go out of scope!!
  }
  ~RandomGenerator() = default;

  double generate() { return m_uniform_distribution(m_generator); }

private:
  std::mt19937_64 m_generator;
  std::uniform_real_distribution<double> m_uniform_distribution;
};

#endif  // GLAC_TESTS_HELPERS_GENERATOR