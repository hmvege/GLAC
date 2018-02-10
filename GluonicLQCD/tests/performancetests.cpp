#include "performancetests.h"

#include <chrono>
#include "config/parameters.h"

// Exponentiation imports
#include "math/exponentiation/expluscher.h"
#include "math/exponentiation/su3exp.h"
#include "math/exponentiation/taylor2exp.h"
#include "math/exponentiation/taylor4exp.h"
#include "math/exponentiation/taylorexp.h"

// Action imports
#include "actions/wilsongaugeaction.h"
#include "actions/luscheraction.h"

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

PerformanceTests::PerformanceTests()
{
    // Sets pointers to use
    m_SU3Generator                      = new SU3MatrixGenerator;
    // Initializing the Mersenne-Twister19937 RNG for the Metropolis algorithm
    m_generator                         = std::mt19937_64(-1);
    m_uniform_distribution              = std::uniform_real_distribution<double>(0,1);

    // Initialize basics, perhaps subclasses for different types of testing?

    // Include 2 derivative speed tests


}

void PerformanceTests::run()
{
    /*
     * Main function for running performance tests.
     */
    m_NTaylorDegree = Parameters::getTaylorPolDegree();
    if (Parallel::Communicator::getProcessRank() == 0) {
        testExponentiationTime(Parameters::getNExpTests());
        testExponentiationAccuracy();
        testRandomGenerators(Parameters::getNRandTests());
    }
    testDerivativeTimeAndAccuracy(Parameters::getNDerivativeTests());
}

void PerformanceTests::testExponentiationTime(int NTests)
{
    printf("\nRunning performance tests for exponentiation timing with %d samples.", NTests);

    // Timers
    double luscherTimer = 0, morningstarTimer = 0, taylor2Timer = 0, taylor4Timer = 0,taylorTimer = 0;
    steady_clock::time_point preUpdate;

    // Initiating different exponentiation matrices
    ExpLuscher expLuscher;
    SU3Exp expMorningstar;
    Taylor2Exp expTaylor2;
    Taylor4Exp expTaylor4;
    TaylorExp expTaylor(m_NTaylorDegree);

    // Generates matrix to test for
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

    // Taylor 16 (numerical precision limit)
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = expTaylor.exp(sampleMatrix);
    }
    taylorTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    printf("\nLuscher exponentiation time:      %8.2f seconds",luscherTimer);
    printf("\nMorningstar exponentiation time:  %8.2f seconds",morningstarTimer);
    printf("\nTaylor2 exponentiation time:      %8.2f seconds",taylor2Timer);
    printf("\nTaylor4 exponentiation time:      %8.2f seconds",taylor4Timer);
    printf("\nTaylor%2d exponentiation time:     %8.2f seconds",m_NTaylorDegree, taylorTimer);
    printf("\nMorningstar/Taylor%2d: %.4f", m_NTaylorDegree, morningstarTimer/taylorTimer);
    printf("\n");
}

void PerformanceTests::testExponentiationAccuracy()
{
    // Initiating different exponentiation matrices
    ExpLuscher expLuscher;
    SU3Exp expMorningstar;
    Taylor2Exp expTaylor2;
    Taylor4Exp expTaylor4;
    TaylorExp expTaylor(m_NTaylorDegree);

    // Generates matrix to test for
    complex temp_diag(0,0);
    SU3 sampleMatrix, temp_sampleMatrix, resultMatrix;
    sampleMatrix = m_SU3Generator->generateRST();
    resultMatrix.zeros();

    // Results to compare for
    double LuscherResult, MorningstarResult, Taylor2Result, Taylor4Result, TaylorNResult;

    // Makes matrix anti hermitian
    temp_sampleMatrix = sampleMatrix.inv();
    temp_sampleMatrix -= sampleMatrix;
    temp_diag = temp_sampleMatrix.trace()/3.0;
    temp_diag.setRe(0);
    temp_sampleMatrix -= temp_diag;
    sampleMatrix = temp_sampleMatrix*0.5;

    // Luscher method
    resultMatrix = expLuscher.exp(sampleMatrix);
    LuscherResult = resultMatrix[0];

    // Morningstar method
    resultMatrix = expMorningstar.exp(sampleMatrix);
    MorningstarResult = resultMatrix[0];

    // Taylor2 method
    resultMatrix = expTaylor2.exp(sampleMatrix);
    Taylor2Result = resultMatrix[0];

    // Taylor4 method
    resultMatrix = expTaylor4.exp(sampleMatrix);
    Taylor4Result = resultMatrix[0];

    // Taylor N method
    resultMatrix = expTaylor.exp(sampleMatrix);
    TaylorNResult = resultMatrix[0];

    // Printing results
    printf("\nComparing first element of test matrix");
    printf("\nLuscher:       %.16f   abs(Luscher - Taylor%2d)     = %.18e", LuscherResult, m_NTaylorDegree, fabs(LuscherResult - TaylorNResult));
    printf("\nMorningstar:   %.16f   abs(Morningstar - Taylor%2d) = %.18e", MorningstarResult, m_NTaylorDegree, fabs(MorningstarResult - TaylorNResult));
    printf("\nTaylor2:       %.16f   abs(Taylor2 - Taylor%2d)     = %.18e", Taylor2Result, m_NTaylorDegree, fabs(Taylor2Result - TaylorNResult));
    printf("\nTaylor4:       %.16f   abs(Taylor4 - Taylor%2d)     = %.18e", Taylor4Result, m_NTaylorDegree, fabs(Taylor4Result - TaylorNResult));
    printf("\nTaylor%2d:      %.16f", m_NTaylorDegree, TaylorNResult);
}

void PerformanceTests::testRandomGenerators(int NTests)
{
    /*
     * Performance tester for the random SU3 matrix generators, RST and random.
     */
    printf("\n\nRunning performance tests for random SU3 matrix generation timing with %d samples.", NTests);

    SU3 RSTRand, FullRand, resultMatrix;

    // Timers
    double RSTTimer = 0, FullRandomTimer = 0;
    steady_clock::time_point preUpdate;

    // RST random
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = RSTRand = m_SU3Generator->generateRandom();
    }
    RSTTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Full random
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        resultMatrix = FullRand = m_SU3Generator->generateRST();
    }
    FullRandomTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    printf("\nRST random generation time:   %8.2f seconds (%8.2E seconds per test)",RSTTimer,RSTTimer/double(NTests));
    printf("\nFull random generation time:  %8.2f seconds (%8.2E seconds per test)",FullRandomTimer,FullRandomTimer/double(NTests));
}

void PerformanceTests::testDerivativeTimeAndAccuracy(int NTests)
{
    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\n\nRunning timing of SU3 derivation methods with %d full lattice derivation tests.",NTests);
    }

    // Timers
    double luscherTimer = 0, morningstarTimer = 0;
    steady_clock::time_point preUpdate;

    // Sets up a full set of sub-lattices
    std::vector<unsigned int> NLatticeDims = {8,8,8,8};
    Parallel::Communicator::setN(NLatticeDims);
    Parallel::Communicator::initializeSubLattice();
    unsigned int subLatticeSize = Parameters::getSubLatticeSize();

    // Initiates the different types of action we will test
    LuscherAction LusAct;
    WilsonGaugeAction MorAct;

    // Generates test lattices and allocates memeory
    Lattice<SU3> *testLattice = new Lattice<SU3>[4];
    Lattice<SU3> *LuscherLattice = new Lattice<SU3>[4];
    Lattice<SU3> *MorningstarLattice = new Lattice<SU3>[4];
    for (int mu = 0; mu < 4; mu++) {
        testLattice[mu].allocate(NLatticeDims);
        LuscherLattice[mu].allocate(NLatticeDims);
        MorningstarLattice[mu].allocate(NLatticeDims);
    }

    // Populates lattices
    for (int mu = 0; mu < 4; mu++) {
        for (unsigned int iSite = 0; iSite < subLatticeSize; iSite++) {
            testLattice[mu][iSite] = m_SU3Generator->generateRST();
        }
        LuscherLattice[mu].zeros();
        MorningstarLattice[mu].zeros();
    }

    // Runs times
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        for (int mu = 0; mu < 4; mu++) {
            MorningstarLattice[mu] = MorAct.getActionDerivative(testLattice,mu);
        }
    }
    morningstarTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Runs timers on the luscher deriative method
    preUpdate = steady_clock::now();
    for (int i = 0; i < NTests; i++) {
        for (int mu = 0; mu < 4; mu++) {
            LuscherLattice[mu] = LusAct.getActionDerivative(testLattice,mu);
        }
    }
    luscherTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Prints results for time test
    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\nLuscher derivation time:      %8.2f seconds (%8.2E seconds per lattice derivative)",luscherTimer,luscherTimer/double(4*NTests));
        printf("\nMorningstar derivation time:  %8.2f seconds (%8.2E seconds per lattice derivative)",morningstarTimer,morningstarTimer/double(4*NTests));
        printf("\nMorningstar/Luscher: %.4f",morningstarTimer/luscherTimer);
        printf("\n");
    }

    // Prints first element of lattice to compare derivative results with
    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\nComparison of the first element of the matrix in the Lattice to each method:");
        printf("\nLuscher:     %.17e",LuscherLattice[0][0][6]);
        printf("\nMorningstar: %.17e",MorningstarLattice[0][0][6]);
        printf("\nAbsolute difference: %.16e",fabs(LuscherLattice[0][0][6] - MorningstarLattice[0][0][6]));
        printf("\n");
    }

    // De-allocates memory
    delete [] testLattice;
    delete [] LuscherLattice;
    delete [] MorningstarLattice;
}
