#ifndef OBSERVABLESIO_H
#define OBSERVABLESIO_H

#include <string>
#include <vector>

namespace IO {
        void writeObservablesToFile(double acceptanceRate,
                                    double averagedObservable,
                                    double varianceObservable,
                                    double stdObservable,
                                    double *observables,
                                    unsigned int NObs,
                                    std::string observableName);
        void writeFlowObservableToFile(double *observables,
                                       std::string observableName,
                                       unsigned int configNumber);
        void writeMatrixToFile(std::vector<double> observables, std::string observableName, unsigned int configNumber, unsigned int N);
}

#endif // OBSERVABLESIO_H
