#ifndef OBSERVABLESIO_H
#define OBSERVABLESIO_H

#include "parameters/parameters.h"
#include "parallelization/communicator.h"

namespace IO {
    class ObservablesIO
    {
    public:
        ObservablesIO();
        void writeObservablesToFile();
        void writeFlowToFile();
    };
}

#endif // OBSERVABLESIO_H
