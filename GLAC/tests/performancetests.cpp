#include "performancetests.h"

#include <chrono>
#include <fstream>
#include <iomanip>

#include "config/parameters.h"
#include "parallelization/communicator.h"

// Exponentiation imports
#include "math/exponentiation/expluscher.h"
#include "math/exponentiation/su3exp.h"
#include "math/exponentiation/taylor2exp.h"
#include "math/exponentiation/taylor4exp.h"
#include "math/exponentiation/taylorexp.h"

// Action imports
#include "actions/wilsongaugeaction.h"
#include "actions/wilsonexplicitder.h"

// Flow import
#include "flow/flow.h"


using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

/*!
 * \brief PerformanceTests::PerformanceTests constructor
 *
 * Initializes random generators and a SU3MatrixGenerator instance.
 */
PerformanceTests::PerformanceTests()
{
    // Sets pointers to use
    m_SU3Generator                      = new SU3MatrixGenerator;
    // Initializing the Mersenne-Twister19937 RNG for the Metropolis algorithm
    m_generator                         = std::mt19937_64(1234);
    m_uniform_distribution              = std::uniform_real_distribution<double>(0,1);
}

/*!
 * \brief PerformanceTests::~PerformanceTests
 *
 * De-allocates the SU3Generator.
 */
PerformanceTests::~PerformanceTests()
{
    delete m_SU3Generator;
}

/*!
 * \brief PerformanceTests::run initiates performance testing on:
 *
 * testExponentiationTime which tests the time of performing matrix exponentiations.
 * testExponentiationAccuracy tests the accuracy of the exponentiation methods.
 * testRandomGenerators tests the random SU3 matrix generators
 * testMatrixMultiplication times the SU3 matrix multiplication method.
 * testDerivativeTimeAndAccuracy tests the lattice derivative methods
 */
void PerformanceTests::run()
{
    /*
     * Main function for running performance tests.
     */
//    m_NTaylorDegree = Parameters::getTaylorPolDegree();
//    if (Parallel::Communicator::getProcessRank() == 0) {
//        testExponentiationTime(Parameters::getNExpTests());
//        testExponentiationAccuracy();
//        testRandomGenerators(Parameters::getNRandTests());
//        testMatrixMultiplication();
//    }
    Parallel::Communicator::setBarrierActive();
//    testDerivativeTimeAndAccuracy(Parameters::getNDerivativeTests());
    testShift();
    testFlow();
}

void PerformanceTests::testExponentiationTime(unsigned int NTests)
{
    /*
     * Method for comparing the timing of different matrix exponentiation methods.
     */
    printf("\nRunning performance tests for exponentiation timing with %d samples.", NTests);

    // Timers
    double luscherTimer = 0, morningstarTimer = 0, taylor2Timer = 0, taylor4Timer = 0;
    steady_clock::time_point preUpdate;

    // Initiating different exponentiation matrices
    ExpLuscher expLuscher;
    SU3Exp expMorningstar;
    Taylor2Exp expTaylor2;
    Taylor4Exp expTaylor4;
    TaylorExp expTaylor(16);

    // Generates matrix to test for
    SU3 sampleMatrix, resultMatrix;
    sampleMatrix = m_SU3Generator->generateRandom();
    resultMatrix.zeros();

    // Luscher
    preUpdate = steady_clock::now();
    for (unsigned int i = 0; i < NTests; i++) {
        resultMatrix = expLuscher.exp(sampleMatrix);
    }
    luscherTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Morningstar
    preUpdate = steady_clock::now();
    for (unsigned int i = 0; i < NTests; i++) {
        resultMatrix = expMorningstar.exp(sampleMatrix);
    }
    morningstarTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Taylor 2
    preUpdate = steady_clock::now();
    for (unsigned int i = 0; i < NTests; i++) {
        resultMatrix = expTaylor2.exp(sampleMatrix);
    }
    taylor2Timer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Taylor 4
    preUpdate = steady_clock::now();
    for (unsigned int i = 0; i < NTests; i++) {
        resultMatrix = expTaylor4.exp(sampleMatrix);
    }
    taylor4Timer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Taylor 16 (numerical precision limit)
    std::vector<double>taylorNtimerResults(m_NTaylorDegree+1);
    for (unsigned int n = 0; n < m_NTaylorDegree+1; n++)
    {
        expTaylor.setTaylorDegree(n);
        preUpdate = steady_clock::now();
        for (unsigned int i = 0; i < NTests; i++) {
            resultMatrix = expTaylor.exp(sampleMatrix);
        }
        taylorNtimerResults[n] = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();
    }

    // Prints timing results
    printf("\nLuscher exponentiation time:      %8.2f seconds",luscherTimer);
    printf("\nMorningstar exponentiation time:  %8.2f seconds",morningstarTimer);
    printf("\nTaylor2 exponentiation time:      %8.2f seconds",taylor2Timer);
    printf("\nTaylor4 exponentiation time:      %8.2f seconds",taylor4Timer);
    printf("\nn TaylorNExpTime Morningstar/TaylorN");
    for (unsigned int n = 0; n < m_NTaylorDegree+1; n++) {
        printf("\n%2d %8.2f %8.4f",n, taylorNtimerResults[n], morningstarTimer/taylorNtimerResults[n]);
    }
    printf("\n");

    // Writes timing results to file
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(16);
        std::ofstream file;
        std::string fname = Parameters::getFilePath()
                          + Parameters::getOutputFolder()
                          + Parameters::getBatchName() + "/"
                          + "observables/" + "exp_timing" + ".dat";
        file.open(fname);
        file << std::fixed << std::setprecision(16);
        file << "Luscher " << luscherTimer << endl;
        file << "Morningstar " << morningstarTimer << endl;
        file << "taylor2 " << taylor2Timer << endl;
        file << "taylor4 " << taylor4Timer << endl;
        file << "n " << "taylorN " << "taylorN/Morningstar" << endl;
        for (unsigned int n = 0; n < m_NTaylorDegree+1; n++) {
            file << n << " " << taylorNtimerResults[n] << " " << morningstarTimer/taylorNtimerResults[n] << endl;
        }
        file.close();
        printf("\n%s written.",fname.c_str());
        std::setprecision(int(oldPrecision));
    }
}

void PerformanceTests::testExponentiationAccuracy()
{
    /*
     * Method for testing the exponentiation accuracy.
     */
    // Initiating different exponentiation matrices
    ExpLuscher expLuscher;
    SU3Exp expMorningstar;
    Taylor2Exp expTaylor2;
    Taylor4Exp expTaylor4;
    TaylorExp expTaylor(m_NTaylorDegree);

    // Generates matrix to test for
    complex temp_diag(0,0);
    SU3 temp_sampleMatrix1, temp_sampleMatrix2, resultMatrix;
    temp_sampleMatrix1 = m_SU3Generator->generateRST();
    resultMatrix.zeros();

    // Results to compare for
    SU3 LuscherResult;
    SU3 MorningstarResult;
    SU3 Taylor2Result;
    SU3 Taylor4Result;
    std::vector<SU3> TaylorNResults(m_NTaylorDegree+1);

    // Makes matrix anti hermitian
    temp_sampleMatrix2 = temp_sampleMatrix1.inv();
    temp_sampleMatrix2 -= temp_sampleMatrix1;
    temp_diag = temp_sampleMatrix2.trace()/3.0;
    temp_diag.setRe(0);
    temp_sampleMatrix2 -= temp_diag;
    const SU3 sampleMatrix = temp_sampleMatrix2*0.5;

    // Luscher method
    resultMatrix = expLuscher.exp(sampleMatrix);
    LuscherResult = resultMatrix;

    // Morningstar method
    resultMatrix = expMorningstar.exp(sampleMatrix);
    MorningstarResult = resultMatrix;

    // Taylor2 method
    resultMatrix = expTaylor2.exp(sampleMatrix);
    Taylor2Result = resultMatrix;

    // Taylor4 method
    resultMatrix = expTaylor4.exp(sampleMatrix);
    Taylor4Result = resultMatrix;

    // Taylor N method
    for (unsigned int n = 0; n < m_NTaylorDegree+1; n++) {
        expTaylor.setTaylorDegree(n);
        resultMatrix = expTaylor.exp(sampleMatrix);
        TaylorNResults[n] = resultMatrix;
    }

    // Element to compare for in matrix
    int cmp_element_index = 0;

    // Comparing results
    printf("\nComparing first element of test matrix exponentiations(absolute error)");
    for (unsigned int n = 0; n < m_NTaylorDegree+1; n++) {
        printf("\n");
        printf("\nLuscher:        %.16f   abs(Luscher - Taylor%2d)     = %.18e", LuscherResult.norm(cmp_element_index), n, fabs(LuscherResult.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)));
        printf("\nMorningstar:    %.16f   abs(Morningstar - Taylor%2d) = %.18e", MorningstarResult.norm(cmp_element_index), n, fabs(MorningstarResult.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)));
        printf("\nMorningstar:    %.16f   rel(Morningstar - Taylor%2d) = %.18e", MorningstarResult.norm(cmp_element_index), n, fabs(MorningstarResult.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index))/fabs(TaylorNResults[n].norm(cmp_element_index)));
        printf("\nTaylor2:        %.16f   abs(Taylor2 - Taylor%2d)     = %.18e", Taylor2Result.norm(cmp_element_index), n, fabs(Taylor2Result.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)));
        printf("\nTaylor4:        %.16f   abs(Taylor4 - Taylor%2d)     = %.18e", Taylor4Result.norm(cmp_element_index), n, fabs(Taylor4Result.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)));
        printf("\nTaylor%2d:      %.16f", n, TaylorNResults[n].norm(cmp_element_index));
        printf("\nTaylor%2d(abs): %.16e", n, fabs(TaylorNResults[n].norm(cmp_element_index) - TaylorNResults[m_NTaylorDegree].norm(cmp_element_index)));
        printf("\nTaylor%2d(rel): %.16e", n, fabs(TaylorNResults[n].norm(cmp_element_index) - TaylorNResults[m_NTaylorDegree].norm(cmp_element_index))/fabs(TaylorNResults[m_NTaylorDegree].norm(cmp_element_index)));
    }

    // Writes timing results to file
    if (Parallel::Communicator::getProcessRank() == 0) {
        auto oldPrecision = cout.precision(16);
        std::ofstream file;
        std::string fname = Parameters::getFilePath()
                          + Parameters::getOutputFolder()
                          + Parameters::getBatchName() + "/"
                          + "observables/" + "exp_precision" + ".dat";
        file.open(fname);
        file << std::fixed << std::setprecision(16);
        file << "luscher " << LuscherResult.norm(cmp_element_index) << endl;
        file << "morningstar " << MorningstarResult.norm(cmp_element_index) << endl;
        file << "taylor2 " << Taylor2Result.norm(cmp_element_index) << endl;
        file << "taylor4 " << Taylor4Result.norm(cmp_element_index) << endl;
        file << "n taylorN "
             << "abs(luscher-taylorN) "
             << "abs(Morningstar-taylorN) "
             << "rel(Morningstar-taylorN) "
             << "abs(Taylor2-taylorN) "
             << "abs(Taylor4-taylorN) "
             << "abs(taylorN-taylor16) "
             << "rel(taylorN-taylor16)" << endl;
        for (unsigned int n = 0; n < m_NTaylorDegree+1; n++) {
            file << n << " ";
            file << TaylorNResults[n].norm(cmp_element_index) << " ";
            file << fabs(LuscherResult.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)) << " ";
            file << fabs(MorningstarResult.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)) << " ";
            file << fabs(MorningstarResult.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index))/fabs(TaylorNResults[n].norm(cmp_element_index)) << " ";
            file << fabs(Taylor2Result.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)) << " ";
            file << fabs(Taylor4Result.norm(cmp_element_index) - TaylorNResults[n].norm(cmp_element_index)) << " ";
            file << fabs(TaylorNResults[n].norm(cmp_element_index) - TaylorNResults[m_NTaylorDegree].norm(cmp_element_index)) << " ";
            file << fabs(TaylorNResults[n].norm(cmp_element_index) - TaylorNResults[m_NTaylorDegree].norm(cmp_element_index))/fabs(TaylorNResults[m_NTaylorDegree].norm(cmp_element_index)) << endl;
        }
        file.close();
        printf("\n%s written.",fname.c_str());
        std::setprecision(int(oldPrecision));
    }
}

void PerformanceTests::testRandomGenerators(unsigned int NTests)
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
    for (unsigned long int i = 0; i < NTests; i++) {
        resultMatrix = RSTRand = m_SU3Generator->generateRandom();
    }
    RSTTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Full random
    preUpdate = steady_clock::now();
    for (unsigned long int i = 0; i < NTests; i++) {
        resultMatrix = FullRand = m_SU3Generator->generateRST();
    }
    FullRandomTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    printf("\nRST random generation time:   %8.2f seconds (%8.2E seconds per test)",RSTTimer,RSTTimer/double(NTests));
    printf("\nFull random generation time:  %8.2f seconds (%8.2E seconds per test)",FullRandomTimer,FullRandomTimer/double(NTests));
}

void PerformanceTests::testDerivativeTimeAndAccuracy(unsigned int NTests)
{
    /*
     * Runs performance tests on the Action derivative times and its accuracy.
     */
    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\n\nRunning timing of SU3 derivation methods with %d full lattice derivation tests.",NTests);
    }


    // Timers
    double luscherTimer = 0, morningstarTimer = 0;
    steady_clock::time_point preUpdate;

    // Sets up a full set of sub-lattices
    std::vector<unsigned int> NLatticeDims = Parameters::getN();
    Parallel::Communicator::initializeSubLattice();
    unsigned long int subLatticeSize = Parameters::getSubLatticeSize();

    // Initiates the different types of action we will test
    WilsonExplicitDer WilsonExpDerAct;
    WilsonGaugeAction WilsonGaugeAct;

    // Generates test lattices and allocates memeory
    Lattice<SU3> *testLattice = new Lattice<SU3>[4];
    Lattice<SU3> *ExplicitExpLattice = new Lattice<SU3>[4];
    Lattice<SU3> *MorningstarLattice = new Lattice<SU3>[4];
    for (int mu = 0; mu < 4; mu++) {
        testLattice[mu].allocate(NLatticeDims);
        ExplicitExpLattice[mu].allocate(NLatticeDims);
        MorningstarLattice[mu].allocate(NLatticeDims);
    }

    // Populates lattices
    for (int mu = 0; mu < 4; mu++) {
        for (unsigned long int iSite = 0; iSite < subLatticeSize; iSite++) {
            testLattice[mu][iSite] = m_SU3Generator->generateRST();
        }
        ExplicitExpLattice[mu].zeros();
        MorningstarLattice[mu].zeros();
    }

    // Runs times
    preUpdate = steady_clock::now();
    for (unsigned long int i = 0; i < NTests; i++) {
        for (int mu = 0; mu < 4; mu++) {
            MorningstarLattice[mu] = WilsonGaugeAct.getActionDerivative(testLattice,mu);
        }
    }
    morningstarTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Runs timers on the luscher deriative method
    preUpdate = steady_clock::now();
    for (unsigned long int i = 0; i < NTests; i++) {
        for (int mu = 0; mu < 4; mu++) {
            ExplicitExpLattice[mu] = WilsonExpDerAct.getActionDerivative(testLattice,mu);
        }
    }
    luscherTimer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    // Prints results for time test
    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\nExplicitExp derivation time:  %8.2f seconds (%8.2E seconds per lattice derivative)",luscherTimer,luscherTimer/double(4*NTests));
        printf("\nMorningstar derivation time:  %8.2f seconds (%8.2E seconds per lattice derivative)",morningstarTimer,morningstarTimer/double(4*NTests));
        printf("\nMorningstar/Luscher: %.4f",morningstarTimer/luscherTimer);
        printf("\n");
    }

    // Prints first element of lattice to compare derivative results with
    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\nComparison of the first element of the matrix in the Lattice to each method:");
        printf("\nExplicitExp: %.17e",ExplicitExpLattice[0][0][6]);
        printf("\nMorningstar: %.17e",MorningstarLattice[0][0][6]);
        printf("\nAbsolute difference: %.16e",fabs(ExplicitExpLattice[0][0][6] - MorningstarLattice[0][0][6]));
        printf("\n");
    }

    // De-allocates memory
    delete [] testLattice;
    delete [] ExplicitExpLattice;
    delete [] MorningstarLattice;
}

void PerformanceTests::testMatrixMultiplication()
{
    /*
     * Runs a SU3 matrix multiplication 10000000 time in order to test the performance.
     */
    unsigned int NTests = 10000000;

    printf("\n\nRunning timing of SU3 matrix multiplication for %d tests",NTests);

    SU3 V0, V1, V2;

    // Timers
    double timer = 0;
    steady_clock::time_point preUpdate;

    preUpdate = steady_clock::now();

    for (unsigned int i = 0; i < NTests; ++i) {
        V0 = m_SU3Generator->generateRandom();
        V1 = m_SU3Generator->generateRandom();
        V2 = V1*V0;
    }

    timer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();

    printf("\nMatrix multiplication tests:  %8.4f seconds (%8.4E seconds per test)",timer,timer/double(NTests));
}

void PerformanceTests::testFlow()
{
    /*
     * Runs performance tests on the flow method.
     */
    Lattice<SU3> *L = new Lattice<SU3>[4];
    for (unsigned int mu = 0; mu < 4; ++mu) {
        L[mu].allocate(m_dim);
    }
    for (unsigned int mu = 0; mu < 4; ++mu) {
        for (unsigned int isite = 0; isite < L[0].m_latticeSize; ++isite) {
            L[mu][isite] = m_SU3Generator->generateRST();
        }
    }

    if (m_processRank == 0) {
        printf("\nRunning flow performance tests for a lattice of size: %d^3 x %d", m_N, m_NT);
    }

    unsigned int NFlow = 100;

    // Timers
    double timerWilsonGauge = 0;
    double timerWilsonExplicit = 0;
    steady_clock::time_point preUpdate;

    // Scopes the test to save memory in the flow method.
    {
        // Runs performance test with Wilson Gauge Action
        WilsonGaugeAction S1;
        Flow flow(&S1);


        preUpdate = steady_clock::now();

        for (unsigned int iflow = 0; iflow < NFlow; ++iflow) {
            flow.flowField(L);
        }

        timerWilsonGauge = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();
        if (m_processRank == 0) {
            printf("\nFlow test for Wilson Action:  %8.4f seconds (%8.4E seconds per flow step)",timerWilsonGauge,timerWilsonGauge/double(NFlow));
        }

    }

    // Resets the lattice
    for (unsigned int mu = 0; mu < 4; ++mu) {
        for (unsigned int isite = 0; isite < L[0].m_latticeSize; ++isite) {
            L[mu][isite] = m_SU3Generator->generateRST();
        }
    }

    // Scopes the test to save memory in the flow method.
    {
        // Runs performance test with Wilson Explicit Exponential Action
        // Sets new flow
        WilsonExplicitDer S2;
        Flow flow2(&S2);

        preUpdate = steady_clock::now();
        for (unsigned int iflow = 0; iflow < NFlow; ++iflow) {
            flow2.flowField(L);
        }

        timerWilsonExplicit = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();
        if (m_processRank == 0) {
            printf("\nFlow test for Wilson explicit Action:  %8.4f seconds (%8.4E seconds per flow step)",timerWilsonExplicit,timerWilsonExplicit/double(NFlow));
        }
    }

    if (m_processRank == 0) {
        printf("\nFlow test improvement for (Wilson Gauge Action) / (Wilson explicit Action) = %8.4f\n", timerWilsonGauge/timerWilsonExplicit);
    }

    delete [] L;
}

void PerformanceTests::testShift()
{
    /*
     * Runs performance tests on the shift method.
     */

    Lattice<SU3> *L = new Lattice<SU3>[4];
    for (unsigned int mu = 0; mu < 4; ++mu) {
        L[mu].allocate(m_dim);
    }
    for (unsigned int mu = 0; mu < 4; ++mu) {
        for (unsigned int isite = 0; isite < L[0].m_latticeSize; ++isite) {
            L[mu][isite] = m_SU3Generator->generateRST();
        }
    }

    if (m_processRank == 0) {
        printf("\nRunning shift performance tests for a lattice of size: %d^3 x %d", m_N, m_NT);
    }

    for (unsigned int mu = 0; mu < 4; ++mu) {
        for (unsigned int isite = 0; isite < L[0].m_latticeSize; ++isite) {
            L[mu][isite] = m_SU3Generator->generateRST();
        }
    }

    unsigned int NTests = 1000;

    // Timers
    double timer = 0;
    steady_clock::time_point preUpdate;
    preUpdate = steady_clock::now();

    for (unsigned int itest = 0; itest < NTests; itest++) {
        for (unsigned int mu = 0; mu < 4; mu++) {
            L[mu] = shift(L[mu], FORWARDS, mu);
            L[mu] = shift(L[mu], BACKWARDS, mu);
        }
    }

    timer = duration_cast<duration<double>>(steady_clock::now() - preUpdate).count();
    if (m_processRank == 0) {
        printf("\nShift test:  %8.4f seconds (%8.4E seconds per flow step)",timer,timer/double(NTests));
    }

    delete [] L;
}
