#include "parameters.h"

// Lattice specific constants
int Parameters::m_NSpatial = 0;
int Parameters::m_NTemporal = 0;
int Parameters::m_latticeSize = 1;
unsigned int Parameters::m_N[] = {0};
int Parameters::m_subLatticeSize = 1;

// Physical lattice size dependent values
double Parameters::m_beta = 0;
double Parameters::m_a = 0;
const double Parameters::r0 = 0.5; // Sommer parameter

// Run specific constants
int Parameters::m_NCf = 0;
int Parameters::m_NCor = 0;
int Parameters::m_NTherm = 0;
int Parameters::m_NUpdates = 0;
int Parameters::m_NFlows = 0;

// Number of points we are storing
int Parameters::m_configSamplePoints = 0;
int Parameters::m_flowSamplePoints = 0;

// Flow
double Parameters::m_flowEpsilon = 0.01;

// Variables holding if we are to calculate and store the thermalization variables
bool Parameters::m_storeThermalizationObservables = false;

// System constants
bool Parameters::m_verbose = true;

// For IO handling
std::string Parameters::m_pwd = "";
std::string Parameters::m_batchName= "";
std::string Parameters::m_inputFolder = "/input/";
std::string Parameters::m_outputFolder = "/output/";

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
     */
    double bval = (beta - 6);
    double a = exp(-1.6805 - 1.7139*bval + 0.8155*bval*bval - 0.6667*bval*bval*bval)*0.5;
    return a/r0;
}

void Parameters::getN(unsigned int *N)
{
    /*
     * Fills the array N with the lattice dimensions.
     */
    for (int i = 0; i < 4; i++)
    {
        N[i] = m_N[i];
    }
}

void Parameters::setNSpatial(int NSpatial)
{
    m_NSpatial = NSpatial;
    m_latticeSize *= NSpatial;
    m_latticeSize *= NSpatial;
    m_latticeSize *= NSpatial;
}

void Parameters::setNTemporal(int NTemporal)
{
    m_NTemporal = NTemporal;
    m_latticeSize *= NTemporal;
}
