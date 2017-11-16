#ifndef FIELDIO_H
#define FIELDIO_H

#include <string>
#include <mpi.h>
#include "math/latticemath.h"
#include "config/parameters.h"
#include "parallelization/index.h"
#include "parallelization/communicator.h"

namespace IO {
    class FieldIO
    {
    private:
        // For writing to file
        static const int m_linkDoubles;// = 72;
        static const int m_linkSize;// = m_linkDoubles*sizeof(double);

        static unsigned int m_N[4];
        static inline double reverseDouble(const double inDouble);
    public:
        FieldIO();
        ~FieldIO();
        static void writeFieldToFile(Links * lattice, int configNumber);
        static void loadFieldConfiguration(std::string filename, Links * lattice);
        static void loadChromaFieldConfiguration(std::string filename, Links *lattice);
    };
}

#endif // FIELDIO_H
