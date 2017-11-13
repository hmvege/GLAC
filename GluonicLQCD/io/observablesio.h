#ifndef OBSERVABLESIO_H
#define OBSERVABLESIO_H

#include "parameters/parameters.h"
#include "parallelization/communicator.h"

namespace IO {
//    class ObservablesIO
//    {
//    public:
//        ObservablesIO();
//        ~ObservablesIO();
        void writeObservablesToFile(double acceptanceRate,
                                    double averagedObservable,
                                    double varianceObservable,
                                    double stdObservable,
                                    double *observables,
                                    std::string observableName);
        void writeObservablesToFile(double acceptanceRate,
                                    double averagedObservable,
                                    double varianceObservable,
                                    double stdObservable,
                                    double *observables,
                                    double *thermalizationObservables,
                                    std::string observableName);
        void writeFlowObservableToFile(double *observableStatistics,
                                       double *t,
                                       double *observables,
                                       std::string observableName);
//    };
}

#endif // OBSERVABLESIO_H
