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
    // IO parameters
    static std::string m_pwd;
    static std::string m_batchName;
    static std::string m_inputFolder;
    static std::string m_outputFolder;

    double calculateLatticeSpacing(double beta);
public:
    Parameters();
    ~Parameters();

    // Setters
    void setNSpatial(int NSpatial) { m_NSpatial = NSpatial; }
    void setNTemporal(int NTemporal) { m_NTemporal = NTemporal; }
    void setLatticeSize(int latticeSize) { m_latticeSize = latticeSize; }
    void setSubLatticeSize(int subLatticeSize) { m_subLatticeSize = subLatticeSize; }
    void setFilePath(std::string pwd) { m_pwd = pwd; }
    void setBatchName(std::string batchName) { m_batchName = batchName; }
    void setInputFolder(std::string inputFolder) { m_inputFolder = inputFolder; }
    void setOutputFolder(std::string outputFolder) { m_outputFolder = outputFolder; }
    void setBeta(double beta);

    // Getters
//    int setNSpatial(int NSpatial) { m_NSpatial = NSpatial; }
//    int setNTemporal(int NTemporal) { m_NTemporal = NTemporal; }
//    int setLatticeSize(int latticeSize) { m_latticeSize = latticeSize; }
//    void setSubLatticeSize(int subLatticeSize) { m_subLatticeSize = subLatticeSize; }
//    void setFilePath(std::string pwd) { m_pwd = pwd; }
//    void setBatchName(std::string batchName) { m_batchName = batchName; }
//    void setInputFolder(std::string inputFolder) { m_inputFolder = inputFolder; }
//    void setOutputFolder(std::string outputFolder) { m_outputFolder = outputFolder; }
//    void setBeta(double beta);

};

#endif // PARAMETERS_H
