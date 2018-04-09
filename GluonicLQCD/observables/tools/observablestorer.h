#ifndef OBSERVABLESTORER_H
#define OBSERVABLESTORER_H

#include <string>
#include <vector>

class ObservableStorer
{
private:
    // Observable name
    std::string m_observableName;
    // Bool to store if we are to normalize the data by number of processors
    bool m_normalizeObservableByProcessor = false;

    // Observable data
    unsigned long m_NObs;
    double m_averagedObservable = 0;
    double m_varianceObservable = 0;
    double m_averagedObservableSquared = 0;
    double m_stdObservable = 0;
    std::vector<double> m_observables;
    std::vector<double> m_observablesSquared;
public:
    ObservableStorer(unsigned long int NSize);
    ~ObservableStorer();

    // Accessor for the observable
    double &operator[](unsigned long int iObs) { return m_observables.at(iObs); }
//    double &operator[](unsigned long int iObs) { return m_observables[iObs]; }

    // Runs statistics, perhaps create its own class? But that increases overhead, so maybe not
    void gatherResults();
    void runStatistics();

    // Printers
    void printStatistics();

    // File writers
    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservableToFile(unsigned long int configNumber);

    // Getters
    std::vector<double> getObservableArray() { return m_observables; }
    double getObservable(unsigned long int iObs) { return m_observables[iObs]; }
    std::string getObservableName() { return m_observableName; }

    // Setters
    void setObservableName(std::string observableName) { m_observableName = observableName; }
    void setNormalizeObservableByProcessor(bool norm) { m_normalizeObservableByProcessor = norm; }
    void reset();
};

#endif // OBSERVABLESTORER_H
