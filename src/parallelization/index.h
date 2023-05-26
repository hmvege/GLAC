/*!
 * \class Index
 *
 * \brief Class for contigious memory accessing.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef INDEX_H
#define INDEX_H

#include <vector>

namespace Parallel {
    class Index
    {
    private:
        static std::vector<unsigned int> m_N;
        static std::vector<long long> m_NTot;
    public:
        Index();
        ~Index();

        // Index getter
        static inline unsigned long cubeIndex(unsigned long int i, unsigned long int j, unsigned long int k, unsigned long int Ni, unsigned long int Nj)
        {
            return i + Ni*(j + Nj*k);
        }

        static inline unsigned long getIndex(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l)
        {
            return i + m_N[0]*(j + m_N[1]*(k + m_N[2]*l));
        }

        static inline long long getGlobalIndex(const long long i, const long long j, const long long k, const long long l)
        {
            return i + m_NTot[0]*(j + m_NTot[1]*(k + m_NTot[2]*l)); // column-major
        }

        // Setters
        static void setN(const std::vector<unsigned int> &N);
        static void setNTot(const unsigned int NSpatial, const unsigned int NTemporal);

        // Getter
        static const std::vector<unsigned int> getN() { return m_N; }
    };
}
#endif // INDEX_H
