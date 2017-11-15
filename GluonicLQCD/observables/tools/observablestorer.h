#ifndef OBSERVABLESTORER_H
#define OBSERVABLESTORER_H

#include <string>
#include "io/observablesio.h"
#include "parallelization/communicator.h"

class ObservableStorer
{
private:
    // Observable name
    std::string m_observableName;
    // Bool to store if we are to normalize the data by number of processors
    bool m_normalizeObservableByProcessor = false;
public:
    ObservableStorer(int NSize);
    ~ObservableStorer();

    // Observable data
    int m_NObs;
    double m_averagedObservable = 0;
    double m_varianceObservable = 0;
    double m_stdObservable = 0;
    double * m_observables; // SLOW COMPARED TO STACK? Making it public is the simplest way to make this work...

    // THIS WILL BE MOVED TO CORRELATOR CLASS, MOST LIKELY!
    void runStatistics();
    // Printers
    void printStatistics();
    // File writers
    void writeObservableToFile();
    void writeFlowObservableToFile(int configNumber);
    // Getters
    double getObservable(int position) { return m_observables[position]; }
    // Setters
    void setObservableName(std::string observableName) {
        m_observableName = observableName;
    }
    void setNormalizeObservableByProcessor(bool norm) { m_normalizeObservableByProcessor = norm; }
};

#endif // OBSERVABLESTORER_H
