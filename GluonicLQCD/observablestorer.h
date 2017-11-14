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
    bool m_procNormalize = false;
public:
    ObservableStorer(int NSize, std::__1::string observableName, bool procNormalize);
    ~ObservableStorer();
    void pushObservable(double newObs, int position);
    void runStatistics();
    void writeObservableToFile();
    void writeFlowObservableToFile(double * flowTime, int configNumber);
};

#endif // OBSERVABLESTORER_H
