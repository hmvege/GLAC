#ifndef FIELDIO_H
#define FIELDIO_H

#include <string>
#include <mpi.h>
#include "math/latticemath.h"
#include "parameters/parameters.h"
#include "parallelization/index.h"
#include "parallelization/communicator.h"

namespace IO {
    class FieldIO
    {
    private:
        // For writing to file
        int m_linkDoubles = 72;
        int m_linkSize = m_linkDoubles*sizeof(double);

        unsigned int m_N[4];
        inline double reverseDouble(const double inDouble);
    public:
        FieldIO();
        void writeFieldToFile(Links * lattice, int configNumber);
        void loadFieldConfiguration(std::string filename, Links * lattice);
        void loadChromaFieldConfiguration(std::string filename, Links *lattice);
    };
}

#endif // FIELDIO_H
