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
        static const int m_SU3Doubles;
        static const int m_SU3Size;

        static std::vector<unsigned int> m_N;
        static inline double reverseDouble(const double inDouble);
        static inline bool check_file_existence (const std::string fname);
    public:
        FieldIO();
        ~FieldIO();
        static void init();
        static void writeFieldToFile(Lattice<SU3> *lattice, int configNumber);
        static void loadFieldConfiguration(std::string filename, Lattice<SU3> *lattice);
        static void loadChromaFieldConfiguration(std::string filename, Lattice<SU3> *lattice);
    };
}

#endif // FIELDIO_H
