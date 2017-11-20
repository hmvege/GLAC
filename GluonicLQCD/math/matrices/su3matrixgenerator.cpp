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

SU3 SU3MatrixGenerator::generateRandom()
{
    /*
     * Generatores a random SU3 matrix.
     * Index map:
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
    H.zeros();
    // Populating the matrix
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            H[6*i + 2*j] = m_uniform_distribution(m_generator);
            H[6*i + 2*j + 1] = m_uniform_distribution(m_generator);
        }
    }
    // Normalizing first column
    double columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H[6*i]*H[6*i] + H[6*i+1]*H[6*i+1];
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[6*i] /= columnLength;
        H[6*i+1] /= columnLength;
    }
    // Using Gram-Schmitt to orthogonalize the second column
    double projectionFactor[2];
    projectionFactor[0] = 0;
    projectionFactor[1] = 0;
    for (int i = 0; i < 3; i++)
    {
        projectionFactor[0] += H[6*i+2]*H[6*i] + H[6*i+3]*H[6*i+1];
        projectionFactor[1] += H[6*i+3]*H[6*i] - H[6*i+2]*H[6*i+1];
    }
    for (int i = 0; i < 3; i++)
    {
        H[6*i+2] -= H[6*i]*projectionFactor[0] - H[6*i+1]*projectionFactor[1];
        H[6*i+3] -= H[6*i]*projectionFactor[1] + H[6*i+1]*projectionFactor[0];
    }
    // Normalizing second column
    columnLength = 0;
    for (int i = 0; i < 3; i++)
    {
        columnLength += H[6*i+2]*H[6*i+2] + H[6*i+3]*H[6*i+3];
    }
    columnLength = sqrt(columnLength);
    for (int i = 0; i < 3; i++)
    {
        H[6*i+2] /= columnLength;
        H[6*i+3] /= columnLength;
    }
    // Taking cross product to produce the last column of our matrix
    H[4] = H[8]*H[12] - H[9]*H[13] - H[14]*H[6] + H[15]*H[7];
    H[5] = H[14]*H[7] + H[15]*H[6] - H[8]*H[13] - H[9]*H[12];
    H[10] = H[14]*H[0] - H[15]*H[1] - H[2]*H[12] + H[3]*H[13];
    H[11] = H[2]*H[13] + H[3]*H[12] - H[14]*H[1] - H[15]*H[0];
    H[16] = H[2]*H[6] - H[3]*H[7] - H[8]*H[0] + H[9]*H[1];
    H[17] = H[8]*H[1] + H[9]*H[0] - H[2]*H[7] - H[3]*H[6];

    return H;
}

SU3 SU3MatrixGenerator::generateRST()
{
    /*
     * Generatores a random SU3 matrix close to unity.
     * Index map:
     * r,s,t =
     * 0 1
     * 2 3
     * R,S,T =
     * 0 1 2
     * 3 4 5
     * 6 7 8
     * H =
     * 0 1 2    0  1   2  3    4  5
     * 3 4 5 =  6  7   8  9   10 11
     * 6 7 8   12 13  14 15   16 17
     */
    // Generates SU2 matrices
    r = generateSU2();
    s = generateSU2();
    t = generateSU2();
    if (m_SU2_uniform_distribution(m_generator) < 0) {
        return RSTMatrixMultiplicationInverse(r,s,t);
    } else {
        return RSTMatrixMultiplication(r,s,t);
    }
}
