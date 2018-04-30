#ifndef SUPERSAMPLER_H
#define SUPERSAMPLER_H

#include "math/lattice.h"
#include "observables/observables.h"
#include <map>

class SuperSampler : public Correlator
{
private:
    const std::string m_observableName = "Weinberg operator";
    const std::string m_observableNameCompact = "weinberg";

    int mu, rho, sigma;
    double m_plaqMultiplicationFactor, m_topcMultiplicationFactor, m_energyMultiplicationFactor, m_wMultiplicationFactor;
    double m_topCharge, m_energy, m_plaquette, m_weinberg;
    Lattice <double> m_tempDiag;
    Lattice<SU3> m_fieldTensorG[6];
    Lattice<SU3> m_clov1, m_clov2, m_U2Temp, m_U3Temp, m_temp;

    // Container for the spatial observables
    std::vector<double> m_topctGatherVector;
    std::vector<double> m_wtGatherVector;
    std::vector<double> m_tempEucl;

    // Creates a object that store the observable
    ObservableStorer * m_plaqObservable = nullptr;
    ObservableStorer * m_topcObservable = nullptr;
    ObservableStorer * m_energyObservable = nullptr;
    ObservableStorer * m_topctObservable = nullptr;
    ObservableStorer * m_wObservable = nullptr; // Weinberg
    ObservableStorer * m_wtObservable = nullptr;

    std::map<int, std::map<int, int>> m_indexMap = {
        {0, {{1, 0}, {2, 1}, {3, 2}}},
        {1, {{2, 4}}},
        {3, {{1, 3}}},
        {2, {{3, 5}}},
    };

//    inline int index_mapper(int i, int j)
//    {
//        /*Function for mapping index onto a lattice array of length 6.*/
//        if (i < j) {
//            return m_indexMap[i][j];
////            return m_indexMap.at(i).at(j);
////            return int(i!=0)*3 + j - 1;
//        } else {
//            return m_indexMap[j][i];
////            return m_indexMap.at(j).at(i);
////            return int(i!=0)*3 + j - 1;
//        }
//    }

    inline int next_index(int i)
    {
        /*Function for getting the next index.*/
        return i % 3 + 1;
    }

public:
    SuperSampler(bool flow);
    ~SuperSampler();
    void calculate(Lattice<SU3> * lattice, unsigned int iObs);
    void initializeObservableStorer(bool storeFlowObservable);

    void writeObservableToFile(double acceptanceRatio);
    void writeFlowObservablesToFile(unsigned int configNumber);
    void reset();
    void runStatistics();
    void printHeader();
    void printObservable(unsigned int iObs);
    void printStatistics();
    std::string getObservableName() { return m_observableName; }
    std::vector<double> getObservablesVector(unsigned int iObs);
    void copyObservable(unsigned int iObs, std::vector<double> obs);
};

#endif // SUPERSAMPLER_H
