#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cmath>
#include <vector>
#include "config/sysprint.h"

class Parameters
{
    friend class SysPrint;
private:
    // Total lattice sizes
    static int m_NSpatial;
    static int m_NTemporal;
    static int m_latticeSize;
    // Sub lattice / parallel related variables
    static unsigned int m_N[4];
    static int m_subLatticeSize;
    static int m_processorsPerDimension[4];
    static bool m_subLatticeSizePreset;
    // Beta value constant
    static double m_beta;
    // Lattice spacing
    static double m_a;
    static const double r0;
    // Program information output
    static bool m_verbose;
    // Variable holding if we are to calculate and store the thermalization variables
    static bool m_storeThermalizationObservables;
    // Variable storing gauge configurations
    static bool m_storeConfigurations;
    // Variable storing if we are to start hot or cold
    static bool m_hotStart;
    // Variable storing what kind if initial hot start we are to use
    static bool m_RSTInit;
    // IO parameters
    static std::string m_pwd;
    static std::string m_batchName;
    static std::string m_inputFolder;
    static std::string m_outputFolder;
    // Run specific variables
    static int m_NCf;
    static int m_NCor;
    static int m_NTherm;
    static int m_NUpdates;
    static int m_NFlows;
    static int m_configSamplePoints;
    static int m_flowSamplePoints;
    // Unit testing
    static bool m_unitTesting;
    static bool m_unitTestingVerbose;
    // Data generation related variables
    static double m_SU3Eps;
    static double m_flowEpsilon;
    static double m_metropolisSeed;
    static double m_randomMatrixSeed;

    // Name of samplers
    static std::string m_expFuncName;
    static std::vector<std::string> m_observablesList;
    static std::vector<std::string> m_flowObservablesList;

    static double calculateLatticeSpacing(double beta);
public:
    Parameters();
    ~Parameters();

//    static void printSystemInfo(); // MOVE INTO OWN CLAS AND MAKE FRIEND?

    /////////////////
    //// Setters ////
    /////////////////
    // Lattice related run variables
    static void setNSpatial(int NSpatial);
    static void setNTemporal(int NTemporal);
    static void setBeta(double beta);
    static void setNCf(int NCf) { m_NCf = NCf; }
    static void setNCor(int NCor) { m_NCor = NCor; }
    static void setNTherm(int NTherm) { m_NTherm = NTherm; }
    static void setNFlows(int NFlows) { m_NFlows = NFlows; }
    static void setNUpdates(int NUpdates) { m_NUpdates = NUpdates; }
    // Data storage related variables
    static void setOutputFolder(std::string outputFolder) { m_outputFolder = outputFolder; }
    static void setInputFolder(std::string inputFolder) { m_inputFolder = inputFolder; }
    static void setStoreConfigurations(bool storeConfigurations) { m_storeConfigurations = storeConfigurations; }
    static void setStoreThermalizationObservables(bool storeThermalizationObservables);
    // Human readable output related variables
    static void setVerbose(bool verbose) { m_verbose = verbose; }
    // Setup related variables
    static void setFilePath(std::string pwd) { m_pwd = pwd; }
    static void setBatchName(std::string batchName) { m_batchName = batchName; }
    static void setHotStart(bool hotStart) { m_hotStart = hotStart; }
    static void setRSTInit(bool RSTInit) { m_RSTInit = RSTInit; }
    // Testing related variables
    static void setUnitTesting(bool unitTesting) { m_unitTesting = unitTesting; }
    static void setUnitTestingVerbose(bool unitTestingVerbose) { m_unitTestingVerbose = unitTestingVerbose; }
    // Data generation related variables
    static void setFlowEpsilon(double flowEpsilon) { m_flowEpsilon = flowEpsilon; }
    static void setSU3Eps(double SU3Eps) { m_SU3Eps = SU3Eps; }
    static void setMetropolisSeed(double metropolisSeed) { m_metropolisSeed = metropolisSeed; }
    static void setRandomMatrixSeed(double randomMatrixSeed) { m_randomMatrixSeed = randomMatrixSeed; }
    // Lattice related variables, initiated after config input
    static void setLatticeSize(int latticeSize) { m_latticeSize = latticeSize; }
    // Sub lattice / parallel related variables
    static void setN(unsigned int *N) { for (int i = 0; i < 4; i++) m_N[i] = N[i]; }
    static void setSubLatticePreset(bool subLatticeSizePreset) { m_subLatticeSizePreset = subLatticeSizePreset; }
    static void setSubLatticeSize(int subLatticeSize) { m_subLatticeSize = subLatticeSize; }
    static void setProcessorsPerDimension(int *processorsPerDimension) { for (int i = 0; i < 4; i++) m_processorsPerDimension[i] = processorsPerDimension[i]; }
    // Internal variables for tracking position in m_observable array
    static void setConfigSamplePoints(int NSamplePoints) { m_configSamplePoints = NSamplePoints; }
    static void setFlowSamplePoints(int NSamplePoints) { m_flowSamplePoints = NSamplePoints; }
    // Name of samplers
    static void setExpFuncName(std::string expFuncName) { m_expFuncName = expFuncName;}
    static void setObservableList(std::vector<std::string> observablesList) { m_observablesList = observablesList; }
    static void setFlowObservablesList(std::vector<std::string> flowObservablesList) { m_flowObservablesList = flowObservablesList; }

    /////////////////
    //// Getters ////
    /////////////////
    // Lattice related run variables
    static int getNSpatial() { return m_NSpatial; }
    static int getNTemporal() { return m_NTemporal; }
    static double getBeta() { return m_beta; }
    static int getNCf() { return m_NCf; }
    static int getNCor() { return m_NCor; }
    static int getNTherm() { return m_NTherm; }
    static int getNUpdates() { return m_NUpdates; }
    static int getNFlows() { return m_NFlows; }
    // Data storage related variables
    static std::string getOutputFolder() { return m_outputFolder; }
    static std::string getInputFolder() { return m_inputFolder; }
    static bool getStoreConfigurations() { return m_storeConfigurations; }
    static bool getStoreThermalizationObservables() { return m_storeThermalizationObservables; }
    // Human readable output related variables
    static bool getVerbose() { return m_verbose; }
    // Setup related variables
    static std::string getFilePath() { return m_pwd; }
    static std::string getBatchName() { return m_batchName; }
    static bool getHotStart() { return m_hotStart; }
    static bool getRSTInit() { return m_RSTInit; }
    // Testing related variables
    static bool getUnitTesting() { return m_unitTesting; }
    static bool getUnitTestingVerbose() { return m_unitTestingVerbose; }
    // Data generation related variables
    static double getFlowEpsilon() { return m_flowEpsilon; }
    static double getSU3Eps() { return m_SU3Eps; }
    static double getMetropolisSeed() { return m_metropolisSeed; }
    static double getRandomMatrixSeed() { return m_randomMatrixSeed; }
    // Lattice related variables, initiated after config input
    static double getLatticeSpacing() { return m_a; }
    static int getLatticeSize() { return m_latticeSize; }
    // Sub lattice / parallel related variables
    static void getN(unsigned int *N) { for (int i = 0; i < 4; i++) N[i] = m_N[i]; }
    static bool getSubLatticePreset() { return m_subLatticeSizePreset; }
    static int getSubLatticeSize() { return m_subLatticeSize; }
    static void getProcessorsPerDimension(int *processorsPerDimension) { for (int i = 0; i < 4; i++) m_processorsPerDimension[i] = processorsPerDimension[i]; }
    // Internal variables for tracking position in m_observable array
    static int getConfigSamplePoints() { return m_configSamplePoints; }
    static int getFlowSamplePoints() { return m_flowSamplePoints; }
    // Name of samplers
    static std::string getExpFuncName() { return m_expFuncName; }
    static std::vector<std::string> getObservablesList() { return m_observablesList; }
    static std::vector<std::string> getFlowObservablesList() { return m_flowObservablesList; }

    // Observable
    // Flow Observable
    // Exponentiation function
//    static

};

#endif // PARAMETERS_H
