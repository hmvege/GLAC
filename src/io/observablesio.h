#ifndef OBSERVABLESIO_H
#define OBSERVABLESIO_H

#include <string>
#include <vector>

namespace IO
{
        void writeObservablesToFile(const double acceptanceRate,
                                    const double averagedObservable,
                                    const double varianceObservable,
                                    const double stdObservable,
                                    const std::vector<double> &observables,
                                    const unsigned int NObs,
                                    const std::string &observableName);
        void writeFlowObservableToFile(const std::vector<double> &observables,
                                       const std::string &observableName,
                                       const unsigned int configNumber);
        void writeMatrixToFile(const std::vector<double> &observables,
                               const std::string &observableName,
                               const unsigned int configNumber,
                               const unsigned int N);
}

#endif // OBSERVABLESIO_H
