#ifndef FIELDIO_H
#define FIELDIO_H

#include <string>
#include "math/links.h"
#include "math/lattice.h"

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
        static void init();
        static void writeFieldToFile(Links * lattice, int configNumber);
        static void loadFieldConfiguration(std::string filename, Links * lattice);
        static void loadChromaFieldConfiguration(std::string filename, Lattice<SU3> *lattice);

        static void loadLatticeFieldConfiguration(std::string filename, Lattice<SU3> *lattice);
    };
}

#endif // FIELDIO_H
