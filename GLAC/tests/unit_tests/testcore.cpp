#include "testcore.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include "io/fieldio.h"

TestCore::TestCore()
{
    /*
     * Base class for unit tests containing matrices and rngs for testing.
     */
    // Initializes parallel usage.
    m_processRank = Parallel::Communicator::getProcessRank();
    m_numprocs = Parallel::Communicator::getNumProc();

    // Initiating the Mersenne Twister random number generator
    m_generator = std::mt19937_64(unsigned(std::time(nullptr) + m_processRank));
    m_uniform_distribution = std::uniform_real_distribution<double>(0,1); // Print out max values to ensure we dont go out of scope!!

    // Initiating the SU3 Matrix generator
    m_SU3Generator = new SU3MatrixGenerator;

    m_verbose = Parameters::getUnitTestingVerbose();

    /////////////////////////////
    //// 3x3 COMPLEX MATRICES ///
    /////////////////////////////
    // Setting up matrix U1
    U1.setComplex(complex(1,1),0);
    U1.setComplex(complex(1,2),2);
    U1.setComplex(complex(1,3),4);
    U1.setComplex(complex(2,1),6);
    U1.setComplex(complex(2,2),8);
    U1.setComplex(complex(2,3),10);
    U1.setComplex(complex(3,1),12);
    U1.setComplex(complex(3,2),14);
    U1.setComplex(complex(3,3),16);
    // Setting up matrix U2
    U2.setComplex(complex(4,4),0);
    U2.setComplex(complex(4,5),2);
    U2.setComplex(complex(4,6),4);
    U2.setComplex(complex(5,4),6);
    U2.setComplex(complex(5,5),8);
    U2.setComplex(complex(5,6),10);
    U2.setComplex(complex(6,4),12);
    U2.setComplex(complex(6,5),14);
    U2.setComplex(complex(6,6),16);
    // Addition
    UAdd.setComplex(complex(5,5),0);
    UAdd.setComplex(complex(5,7),2);
    UAdd.setComplex(complex(5,9),4);
    UAdd.setComplex(complex(7,5),6);
    UAdd.setComplex(complex(7,7),8);
    UAdd.setComplex(complex(7,9),10);
    UAdd.setComplex(complex(9,5),12);
    UAdd.setComplex(complex(9,7),14);
    UAdd.setComplex(complex(9,9),16);
    // Subtraction
    USub.setComplex(complex(-3,-3),0);
    USub.setComplex(complex(-3,-3),6);
    USub.setComplex(complex(-3,-3),12);
    USub.setComplex(complex(-3,-3),2);
    USub.setComplex(complex(-3,-3),8);
    USub.setComplex(complex(-3,-3),14);
    USub.setComplex(complex(-3,-3),4);
    USub.setComplex(complex(-3,-3),10);
    USub.setComplex(complex(-3,-3),16);
    // Multiplication
    UMul.setComplex(complex(-9,44),0);
    UMul.setComplex(complex(-15,47),2);
    UMul.setComplex(complex(-21,50),4);
    UMul.setComplex(complex(6,56),6);
    UMul.setComplex(complex(0,62),8);
    UMul.setComplex(complex(-6,68),10);
    UMul.setComplex(complex(21,68),12);
    UMul.setComplex(complex(15,77),14);
    UMul.setComplex(complex(9,86),16);
    // Conjugation
    UC.setComplex(complex(1,-1),0);
    UC.setComplex(complex(1,-2),2);
    UC.setComplex(complex(1,-3),4);
    UC.setComplex(complex(2,-1),6);
    UC.setComplex(complex(2,-2),8);
    UC.setComplex(complex(2,-3),10);
    UC.setComplex(complex(3,-1),12);
    UC.setComplex(complex(3,-2),14);
    UC.setComplex(complex(3,-3),16);
    // Transpose
    UT.setComplex(complex(1,1),0);
    UT.setComplex(complex(2,1),2);
    UT.setComplex(complex(3,1),4);
    UT.setComplex(complex(1,2),6);
    UT.setComplex(complex(2,2),8);
    UT.setComplex(complex(3,2),10);
    UT.setComplex(complex(1,3),12);
    UT.setComplex(complex(2,3),14);
    UT.setComplex(complex(3,3),16);
    // Conjugate transpose
    UCT.setComplex(complex(1,-1),0);
    UCT.setComplex(complex(2,-1),2);
    UCT.setComplex(complex(3,-1),4);
    UCT.setComplex(complex(1,-2),6);
    UCT.setComplex(complex(2,-2),8);
    UCT.setComplex(complex(3,-2),10);
    UCT.setComplex(complex(1,-3),12);
    UCT.setComplex(complex(2,-3),14);
    UCT.setComplex(complex(3,-3),16);
    // RST Matrix multiplier
    U_RST.setComplex(complex(-8,9),0);
    U_RST.setComplex(complex(-123,82),2);
    U_RST.setComplex(complex(-78,83),4);
    U_RST.setComplex(complex(-20,37),6);
    U_RST.setComplex(complex(-391,342),8);
    U_RST.setComplex(complex(-234,319),10);
    U_RST.setComplex(complex(2,6),12);
    U_RST.setComplex(complex(48,58),14);
    U_RST.setComplex(complex(44,35),16);
    // r SU2 matrix in RST multiplier test
    s_r.setComplex(complex(1,2),0);
    s_r.setComplex(complex(3,4),2);
    s_r.setComplex(complex(5,6),4);
    s_r.setComplex(complex(7,8),6);
    // s SU2 matrix in RST multiplier test
    s_s.setComplex(complex(2,5),0);
    s_s.setComplex(complex(4,8),2);
    s_s.setComplex(complex(2,6),4);
    s_s.setComplex(complex(10,3),6);
    // t SU2 matrix in RST multiplier test
    s_t.setComplex(complex(7,2),0);
    s_t.setComplex(complex(6,1),2);
    s_t.setComplex(complex(6,4),4);
    s_t.setComplex(complex(5,2),6);
    // Traced matrix multiplication
    m_tracedMatrix = -21.0; // re(tr(U1*U3))
    UTrace.setComplex(complex(4,4),0);
    UTrace.setComplex(complex(5,5),2);
    UTrace.setComplex(complex(2,6),4);
    UTrace.setComplex(complex(6,4),6);
    UTrace.setComplex(complex(5,5),8);
    UTrace.setComplex(complex(1,6),10);
    UTrace.setComplex(complex(6,4),12);
    UTrace.setComplex(complex(9,5),14);
    UTrace.setComplex(complex(2,6),16);
    // Hermitian U1 matrix
    U1Hermitian.setComplex(complex(-4,0),0);
    U1Hermitian.setComplex(complex(-3,1),2);
    U1Hermitian.setComplex(complex(-2,2),4);
    U1Hermitian.setComplex(complex(-3,-1),6);
    U1Hermitian.setComplex(complex(-2,0),8);
    U1Hermitian.setComplex(complex(-1,1),10);
    U1Hermitian.setComplex(complex(-2,-2),12);
    U1Hermitian.setComplex(complex(-1,-1),14);
    U1Hermitian.setComplex(complex(0,0),16);
    // Anti-hermitian U1 matrix
    U1AntiHermitian.setComplex(complex(0,-4),0);
    U1AntiHermitian.setComplex(complex(-1,-3),2);
    U1AntiHermitian.setComplex(complex(-2,-2),4);
    U1AntiHermitian.setComplex(complex(1,-3),6);
    U1AntiHermitian.setComplex(complex(0,-2),8);
    U1AntiHermitian.setComplex(complex(-1,-1),10);
    U1AntiHermitian.setComplex(complex(2,-2),12);
    U1AntiHermitian.setComplex(complex(1,-1),14);
    U1AntiHermitian.setComplex(complex(0,0),16);

    /////////////////////////////
    //// 2x2 COMPLEX MATRICE ////
    /////////////////////////////
    // s1
    s1.setComplex(complex(1,1),0);
    s1.setComplex(complex(1,2),2);
    s1.setComplex(complex(2,1),4);
    s1.setComplex(complex(2,2),6);
    // s2
    s2.setComplex(complex(4,4),0);
    s2.setComplex(complex(4,5),2);
    s2.setComplex(complex(5,4),4);
    s2.setComplex(complex(5,5),6);
    // Adding
    sAdd.setComplex(complex(5,5),0);
    sAdd.setComplex(complex(5,7),2);
    sAdd.setComplex(complex(7,5),4);
    sAdd.setComplex(complex(7,7),6);
    // Subtracting
    sSub.setComplex(complex(-3,-3),0);
    sSub.setComplex(complex(-3,-3),2);
    sSub.setComplex(complex(-3,-3),4);
    sSub.setComplex(complex(-3,-3),6);
    // Multiplying
    sMul.setComplex(complex(-3,22),0);
    sMul.setComplex(complex(-6,24),2);
    sMul.setComplex(complex(6,30),4);
    sMul.setComplex(complex(3,34),6);
    // Conjugate of s1
    sC.setComplex(complex(1,-1),0);
    sC.setComplex(complex(1,-2),2);
    sC.setComplex(complex(2,-1),4);
    sC.setComplex(complex(2,-2),6);
    // Transpose of s1
    sT.setComplex(complex(1,1),0);
    sT.setComplex(complex(2,1),2);
    sT.setComplex(complex(1,2),4);
    sT.setComplex(complex(2,2),6);
    // Complex conjugate of s1
    sCT.setComplex(complex(1,-1),0);
    sCT.setComplex(complex(2,-1),2);
    sCT.setComplex(complex(1,-2),4);
    sCT.setComplex(complex(2,-2),6);

    /////////////////////////////
    //// Complex unit testing ///
    /////////////////////////////
    z1 = complex(1,2);
    z2 = complex(3,4);
    // Adding
    zAdd = complex(4,6);
    // Subtracting
    zSub = complex(-2,-2);
    // Multiplying
    zMul = complex(-5,10);
    // Division
    zDiv = complex(0.44,0.08);
    // Conjugate 1, conjugate()/c()
    zConj = complex(1,-2);
    // Norm
    zNorm = 2.23606797749979;
    // Norm squared
    zNormSquared = 5;
    // Set to minus operator
    zSetToMinus = complex(-1,-2);

    /////////////////////////////
    //// Lattice operations /////
    /////////////////////////////
    latticeDoubleValue1 = 2.5;
    latticeDoubleValue2 = 1.5;
    m_N = Parameters::getNSpatial();
    m_NT = Parameters::getNTemporal();
//    Parameters::setSubLatticePreset(false);
    Parallel::Communicator::initializeSubLattice();
    m_subLatticeSize = Parameters::getSubLatticeSize();
    m_dim = Parameters::getN();
    IO::FieldIO::init();
    latticeSU3_U1.allocate(m_dim);
    latticeSU3_U2.allocate(m_dim);
    latticeComplex_z1.allocate(m_dim);
    latticeComplex_z2.allocate(m_dim);
    latticeDouble1.allocate(m_dim);
    latticeDouble2.allocate(m_dim);
    for (unsigned int iSite = 0; iSite < latticeSU3_U1.m_latticeSize; iSite++) {
        latticeSU3_U1[iSite] = U1;
        latticeSU3_U2[iSite] = U2;
        latticeComplex_z1[iSite] = z1;
        latticeComplex_z2[iSite] = z2;
        latticeDouble1[iSite] = latticeDoubleValue1;
        latticeDouble2[iSite] = latticeDoubleValue2;
    }
}

TestCore::~TestCore() {
    delete m_SU3Generator;
}
