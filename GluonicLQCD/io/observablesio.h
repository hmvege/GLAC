#ifndef OBSERVABLESIO_H
#define OBSERVABLESIO_H

#include "config/parameters.h"

namespace IO {
        void writeObservablesToFile(double acceptanceRate,
                                    double averagedObservable,
                                    double varianceObservable,
                                    double stdObservable,
                                    double *observables,
                                    int NObs,
                                    std::string observableName);
        void writeFlowObservableToFile(double *observables,
                                       std::string observableName,
                                       int configNumber);
}

#endif // OBSERVABLESIO_H
