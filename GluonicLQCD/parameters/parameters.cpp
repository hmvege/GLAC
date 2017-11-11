#include "parameters.h"

int Parameters::m_NSpatial = 0;
int Parameters::m_NTemporal = 0;
int Parameters::m_latticeSize = 0;
unsigned int Parameters::m_N[] = {0};
int Parameters::m_subLatticeSize = 0;
double Parameters::m_beta = 0;
double Parameters::m_a = 0;

std::string Parameters::m_pwd = "";
std::string Parameters::m_filename = "";
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
    m_beta = beta;
    m_a = calculateLatticeSpacing(beta;
}

double Parameters::calculateLatticeSpacing(double beta)
{
    double r = 0.5;
    double bval = (beta - 6);
    double a = exp(-1.6805 - 1.7139*bval + 0.8155*bval*bval - 0.6667*bval*bval*bval)*0.5;
    return a;
}
