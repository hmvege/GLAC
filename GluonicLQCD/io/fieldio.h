/*!
 * \class FieldIO
 *
 * \brief Class for reading and writing lattice to file.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
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
        static const long long m_linkDoubles;// = 72;
        static const long long m_linkSize;// = m_linkDoubles*sizeof(double);
        static const long long m_SU3Doubles;
        static const long long m_SU3Size;

        static std::vector<unsigned int> m_N;
        static inline double reverseDouble(const double inDouble);
        static inline bool check_file_existence (const std::string fname);
    public:
        FieldIO();
        ~FieldIO();
        static void init();
        static void writeFieldToFile(Lattice<SU3> *lattice, unsigned int configNumber);
        static void writeDoublesFieldToFile(Lattice<double> lattice, unsigned int configNumber, std::string observable);
        static void loadFieldConfiguration(std::string filename, Lattice<SU3> *lattice);
        static void loadChromaFieldConfiguration(std::string filename, Lattice<SU3> *lattice);
    };
}

#endif // FIELDIO_H
