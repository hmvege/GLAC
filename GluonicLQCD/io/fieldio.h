#ifndef FIELDIO_H
#define FIELDIO_H

#include "math/latticemath.h"

namespace IO {
    class FieldIO
    {
    public:
        FieldIO();
        void writeFieldToFile(Links * lattice);
    };
}

#endif // FIELDIO_H
