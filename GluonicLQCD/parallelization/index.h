#ifndef INDEX_H
#define INDEX_H

#include <vector>

namespace Parallel {
    class Index
    {
    private:
        static std::vector<unsigned int> m_N;
        static unsigned int m_NTot[4];
    public:
        Index();
        ~Index();

        // Index getter
        inline static unsigned int cubeIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int Ni, unsigned int Nj) {
            return i + Ni*(j + Nj*k);
        }

        static unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
            return i + m_N[0]*(j + m_N[1]*(k + m_N[2]*l));
        }

        static unsigned int getGlobalIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
            return i + m_NTot[0]*(j + m_NTot[1]*(k + m_NTot[2]*l)); // column-major
        }

        // Setters
        static void setN(std::vector<unsigned int> N);
        static void setNTot(int NSpatial, int NTemporal);
    };
}
#endif // INDEX_H
