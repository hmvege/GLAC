#include "performancetests.h"

#include <chrono>
#include "config/parameters.h"
#include "math/exponentiation/expluscher.h"
#include "math/exponentiation/su3exp.h"
#include "math/exponentiation/taylor2exp.h"
#include "math/exponentiation/taylor4exp.h"

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

PerformanceTests::PerformanceTests()
{
    // Sets pointers to use
    m_SU3Generator                      = new SU3MatrixGenerator;
    // Initializing the Mersenne-Twister19937 RNG for the Metropolis algorithm
    m_generator                         = std::mt19937_64(Parameters::getMetropolisSeed());
    m_uniform_distribution              = std::uniform_real_distribution<double>(0,1);

    // Initialize basics, perhaps subclasses for different types of testing?

    // Include 4 tests with exponentiation tests

    // Include 2 derivative speed tests

    // Include RST, Random tests

}

void PerformanceTests::testExponentiation(int NTests)
{
    printf("\nRunning performance tests for exponentiation timing.");

    // Timers
    double luscherTimer = 0, morningstarTimer = 0, taylor2Timer = 0, taylor4Timer = 0;
    steady_clock::time_point preUpdate;

    ExpLuscher expLuscher;
    SU3Exp expMorningstar;
    Taylor2Exp expTaylor2;
    Taylor4Exp expTaylor4;

    SU3 sampleMatrix, resultMatrix;
    sampleMatrix = m_SU3Generator->generateRandom();
    resultMatrix.zeros();

    // Luscher
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = expLuscher.exp(sampleMatrix);
    }
    luscherTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Morningstar
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = expMorningstar.exp(sampleMatrix);
    }
    morningstarTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Taylor 2
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = expTaylor2.exp(sampleMatrix);
    }
    taylor2Timer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Taylor 4
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = expTaylor4.exp(sampleMatrix);
    }
    taylor4Timer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    printf("\nLuscher time for %d exponentiation: %.2f",NTests,luscherTimer);
    printf("\nMorningstar time for %d exponentiation: %.2f",NTests,morningstarTimer);
    printf("\nTaylor2 time for %d exponentiation: %.2f",NTests,taylor2Timer);
    printf("\nTaylor4 time for %d exponentiation: %.2f",NTests,taylor4Timer);
}


