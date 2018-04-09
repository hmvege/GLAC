#include "parameters.h"
#include <cmath>

// Lattice specific constants
unsigned int Parameters::m_NSpatial = 0;
unsigned int Parameters::m_NTemporal = 0;
unsigned int Parameters::m_latticeSize = 1;

// Sub lattice / parallel related variables
std::vector<unsigned int> Parameters::m_N;
unsigned int Parameters::m_subLatticeSize = 1;
int Parameters::m_processorsPerDimension[4] = {0,0,0,0};
bool Parameters::m_subLatticeSizePreset = false;

// Physical lattice size dependent values
double Parameters::m_beta = 0;
double Parameters::m_a = 0;
const double Parameters::r0 = 0.5; // Sommer parameter

// Run specific constants
unsigned int Parameters::m_NCf = 0;
unsigned int Parameters::m_NCor = 0;
unsigned int Parameters::m_NTherm = 0;
unsigned int Parameters::m_NUpdates = 0;
unsigned int Parameters::m_NFlows = 0;

// Number of points we are storing
unsigned int Parameters::m_configSamplePoints = 0;
unsigned int Parameters::m_flowSamplePoints = 0;

// Variables holding if we are to calculate and store the thermalization variables
bool Parameters::m_storeThermalizationObservables = false;

// Variable storing gauge configurations
bool Parameters::m_storeConfigurations = false;

// Variable storing if we are to start hot or cold
bool Parameters::m_hotStart = false;

// Variable storing what kind if initial hot start we are to use
bool Parameters::m_RSTHotStart = false;

// System constants
bool Parameters::m_verbose = true;

// For IO handling
std::string Parameters::m_pwd = "";
std::string Parameters::m_batchName= "";
std::string Parameters::m_inputFolder = "/input/";
std::string Parameters::m_outputFolder = "/output/";

// Unit testing
bool Parameters::m_unitTesting = false;
bool Parameters::m_unitTestingVerbose = false;
bool Parameters::m_testLatticeGaugeInvariance = false;
std::string Parameters::m_latticeFileNameToCheck = "";

// Performance testing
bool Parameters::m_performanceTesting = false;
int Parameters::m_NExpTests = 0;
int Parameters::m_NRandTests = 0;
int Parameters::m_NDerivativeTests = 0;
int Parameters::m_NTaylorPolDegree = 8;

// Data generation related variables
double Parameters::m_SU3Eps = 0.24;
double Parameters::m_flowEpsilon = 0.01;
double Parameters::m_metropolisSeed = 0;
double Parameters::m_randomMatrixSeed = 0;

// Name of samplers
std::string Parameters::m_expFuncName = "";
std::vector<std::string> Parameters::m_observablesList = {};
std::vector<std::string> Parameters::m_flowObservablesList = {};

// Field configurations
bool Parameters::m_loadFieldConfigurations = false;
bool Parameters::m_loadChromaConfigurations = false;
std::vector<std::string> Parameters::m_fieldConfigurationFileNames = {};

// Variable for storing if we are to load and run from a given configuration
bool Parameters::m_loadConfigAndRun = false;
int Parameters::m_configStartNumber = 0;

// Parameter for storing the sampling frequency of the field
int Parameters::m_samplingFrequency = 25;

// Debug parameter
bool Parameters::m_debug = false;

Parameters::Parameters()
{

}

Parameters::~Parameters()
{

}

void Parameters::setBeta(double beta)
{
    /*
     * Sets the beta value, and finds the lattice spacing a.
     */
    m_beta = beta;
    m_a = calculateLatticeSpacing(beta);
}

double Parameters::calculateLatticeSpacing(double beta)
{
    /*
     * Calculates the lattice spaing from the beta value.
     * From paper by M. Guagnelli, arXiv:hep-lat/9806005
     */
    double bval = (beta - 6);
    double a = exp(-1.6805 - 1.7139*bval + 0.8155*bval*bval - 0.6667*bval*bval*bval)*0.5;
    return a/r0;
}

void Parameters::setNSpatial(unsigned long int NSpatial)
{
    m_NSpatial = NSpatial;
    m_latticeSize *= (unsigned long int) NSpatial*NSpatial*NSpatial;
}

void Parameters::setNTemporal(unsigned long int NTemporal)
{
    m_NTemporal = NTemporal;
    m_latticeSize *= (unsigned long int) NTemporal;
}

void Parameters::setMetropolisSeed(double metropolisSeed)
{
        metropolisSeed += double(Parallel::Communicator::getNumProc()) + double(Parallel::Communicator::getProcessRank());
        m_metropolisSeed = metropolisSeed;
}

void Parameters::setRandomMatrixSeed(double randomMatrixSeed)
{
    randomMatrixSeed += double(Parallel::Communicator::getProcessRank());
    m_randomMatrixSeed = randomMatrixSeed;
}
