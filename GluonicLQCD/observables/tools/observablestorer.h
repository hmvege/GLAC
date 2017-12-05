#ifndef OBSERVABLESTORER_H
#define OBSERVABLESTORER_H

#include <string>

class ObservableStorer
{
private:
    // Observable name
    std::string m_observableName;
    // Bool to store if we are to normalize the data by number of processors
    bool m_normalizeObservableByProcessor = false;

    // Observable data
    int m_NObs;
    double m_averagedObservable = 0;
    double m_varianceObservable = 0;
    double m_averagedObservableSquared = 0;
    double m_stdObservable = 0;
    double * m_observables; // SLOW COMPARED TO STACK? --> OVERLOAD THIS!!!
    double * m_observablesSquared;
public:
    ObservableStorer(int NSize);
    ~ObservableStorer();

    // Accessor for the observable
    double &operator[](int iObs) { return m_observables[iObs]; }

    // Runs statistics, perhaps create its own class? But that increases overhead, so maybe not
    void gatherResults();
    void runStatistics();
    // Printers
    void printStatistics();

    // File writers
    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservableToFile(int configNumber);
    // Getters
    double getObservable(int position) { return m_observables[position]; }
    std::string getObservableName() { return m_observableName; }
    // Setters
    void setObservableName(std::string observableName) { m_observableName = observableName; }
    void setNormalizeObservableByProcessor(bool norm) { m_normalizeObservableByProcessor = norm; }
    void reset();
};

#endif // OBSERVABLESTORER_H
