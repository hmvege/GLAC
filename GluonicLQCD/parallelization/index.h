#ifndef INDEX_H
#define INDEX_H

namespace Parallel {
    class Index
    {
    private:
        static unsigned int m_N[4];
        static unsigned int m_NTot[4];
    public:
        Index();
        ~Index();

        // Index getter
        static unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l);
        static unsigned int getGlobalIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l);

        // Setters
        static void setN(unsigned int *N);
        static void setNTot(int NSpatial, int NTemporal);
    };
}
#endif // INDEX_H
