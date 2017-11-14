#ifndef OBSERVABLESTORER_H
#define OBSERVABLESTORER_H

#include <string>
#include "io/observablesio.h"
#include "parallelization/communicator.h"

class ObservableStorer
{
private:
    // Observable data
    int m_NSize = 0;
    double * m_observables;
    double * m_observablesSquared;
    double m_averagedObservable = 0;
    double m_varianceObservable = 0;
    double m_stdObservable = 0;
    // Observable name
    std::string m_observableName;
    // Bool to store if we are to normalize the data by number of processors
    bool m_normalizeObservableByProcessor = false;
public:
    ObservableStorer(int m_NSize);
    ~ObservableStorer();
    void pushObservable(double newObs, int position);
    void runStatistics();
    // Printers
    void printStatistics();
    // File writers
    void writeObservableToFile();
    void writeFlowObservableToFile(double * flowTime, int configNumber);
    // Getters
    double getObservable(int position) { return m_observables[position]; }
    // Setters
    void setObservableName(std::string observableName) { m_observableName = observableName; }
    void setNormalizeObservableByProcessor(bool norm) { m_normalizeObservableByProcessor = norm; }
};

#endif // OBSERVABLESTORER_H
