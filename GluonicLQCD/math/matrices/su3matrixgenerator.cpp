#include "su3matrixgenerator.h"
#include <random>
#include <iomanip>
#include <iostream>
#include "math/functions.h"
#include "complex.h"
#include "config/parameters.h"

// For testing
#include "tests/unittests.h"

using std::cout;
using std::endl;

/*
 * Function for generating random SU3 matrices
 */

SU3MatrixGenerator::SU3MatrixGenerator()
{
    m_epsilon = Parameters::getSU3Eps();
    m_epsilonSquared = m_epsilon*m_epsilon;
    m_sqrtOneMinusEpsSquared = sqrt(1 - m_epsilonSquared);
    // Initializes RNGs
    m_generator = std::mt19937_64 (Parameters::getRandomMatrixSeed());
    m_uniform_distribution = std::uniform_real_distribution<double> (-m_epsilon,m_epsilon);
    m_SU2_uniform_distribution = std::uniform_real_distribution<double> (-0.5,0.5);
    // Ensures RST and regular matrices start at zero.
    H.zeros();
    X.zeros();
    // Populates Pauli matrices
    sigma[0].mat[2] = 1;
    sigma[0].mat[4] = 1;
    sigma[1].mat[3] = -1;
    sigma[1].mat[5] = 1;
    sigma[2].mat[0] = 1;
    sigma[2].mat[6] = -1;
}

SU3MatrixGenerator::~SU3MatrixGenerator()
{
}
