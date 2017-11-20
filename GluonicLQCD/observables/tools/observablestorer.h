#ifndef OBSERVABLESTORER_H
#define OBSERVABLESTORER_H

#include <string>
#include "io/observablesio.h"
#include "parallelization/communicator.h"

struct ObservableStorer
{
    ObservableStorer(int NSize);
    ~ObservableStorer();

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
    double * m_observables; // SLOW COMPARED TO STACK?
    double * m_observablesSquared;

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
    // Setters
    void setObservableName(std::string observableName) { m_observableName = observableName; }
    void setNormalizeObservableByProcessor(bool norm) { m_normalizeObservableByProcessor = norm; }
};

#endif // OBSERVABLESTORER_H
