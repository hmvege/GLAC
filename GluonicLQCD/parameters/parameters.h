#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cmath>

class Parameters
{
private:
    // Total lattice sizes
    static int m_NSpatial;
    static int m_NTemporal;
    static int m_latticeSize;
    // Sub lattice sizes
    static unsigned int m_N[4];
    static int m_subLatticeSize;
    // Beta value constant
    static double m_beta;
    // Lattice spacing
    static double m_a;
    static const double r0;
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

    // Flow
    static double m_flowEpsilon; // TEMP SET TO 0.02!

    static double calculateLatticeSpacing(double beta);
public:
    Parameters();
    ~Parameters();

    // Setters
    static void setNSpatial(int NSpatial);
    static void setNTemporal(int NTemporal);
    static void setNCf(int NCf) { m_NCf = NCf; }
    static void setNCor(int NCor) { m_NCor = NCor; }
    static void setNTherm(int NTherm) { m_NTherm = NTherm; }
    static void setNUpdates(int NUpdates) { m_NUpdates = NUpdates; }
    static void setNFlows(int NFlows) { m_NFlows = NFlows; }
    static void setLatticeSize(int latticeSize) { m_latticeSize = latticeSize; }
    static void setSubLatticeSize(int subLatticeSize) { m_subLatticeSize = subLatticeSize; }
    static void setFilePath(std::string pwd) { m_pwd = pwd; }
    static void setBatchName(std::string batchName) { m_batchName = batchName; }
    static void setInputFolder(std::string inputFolder) { m_inputFolder = inputFolder; }
    static void setOutputFolder(std::string outputFolder) { m_outputFolder = outputFolder; }
    static void setBeta(double beta);

    // Getters
    static int getNSpatial() { return m_NSpatial; }
    static int getNTemporal() { return m_NTemporal; }
    static int getLatticeSize() { return m_latticeSize; }
    static int getSubLatticeSize() { return m_subLatticeSize; }
    static int getNCf() { return m_NCf; }
    static int getNCor() { return m_NCor; }
    static int getNTherm() { return m_NTherm; }
    static int getNUpdates() { return m_NUpdates; }
    static int getNFlows() { return m_NFlows; }
    static double getFlowEpsilon() { return m_flowEpsilon; }
    static double getBeta() { return m_beta; }
    static double getLatticeSpacing() { return m_a; }
    static void getN(unsigned int *N);
    static std::string getFilePath() { return m_pwd; }
    static std::string getBatchName() { return m_batchName; }
    static std::string getInputFolder() { return m_inputFolder; }
    static std::string getOutputFolder() { return m_outputFolder; }

};

#endif // PARAMETERS_H
