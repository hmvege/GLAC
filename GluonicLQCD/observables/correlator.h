#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "math/links.h"
#include "parallelization/index.h"
//#include "parallelization/neighbours.h"
#include "parallelization/communicator.h"
#include "tools/observablestorer.h"
#include <vector>
#include <iostream>


using std::cout;
using std::endl;

class Correlator
{
protected:
    // Lattice dimension array
    unsigned int *m_N;
    // (sub)lattice size
    double m_latticeSize;
    // Lattice spacing
    double m_a;
    // Lorentz indices
    int muIndex[4];
    int nuIndex[4];
    // Position vector for handling the shift-method in parallelization
    std::vector<int> m_position;
    // Observable name
    const int m_headerWidth = 20;
    const std::string m_observableName = "Correlator";
    const std::string m_observableNameCompact = "correlator";
    // Creates a object that store the observable
    ObservableStorer * m_observable = nullptr;
    // Object that holds info of whether or not we are storing flow observables
    bool m_storeFlowObservable = false;
    // Function for updating the mu and nu indices
    inline void updateMuIndex(int mu) {
        for (int i = 0; i < 4; i++)
        {
            muIndex[i] = 0;
        }
        muIndex[mu] = 1;
    }
    inline void updateNuIndex(int nu) {
        for (int i = 0; i < 4; i++)
        {
            nuIndex[i] = 0;
        }
        nuIndex[nu] = 1;
    }
    // Function for accessing the clover index
    inline int cloverIndex(int mu, int nu)
    {
        /*
         * Used for accessing the clover SU3 elements in an contigious fashion.
         */
        return (4*mu + nu);
    }
public:
    Correlator(bool storeFlowObservable);
    virtual ~Correlator();
    virtual void calculate(Links *lattice, int iObs);
    virtual void calculate(SU3 *U, int iObs);
    virtual void writeStatisticsToFile(double acceptanceRatio);
    virtual void writeFlowObservablesToFile(int iFlow);
    virtual void runStatistics();

    // Printers
    virtual void printObservable(int iObs);
    virtual void printHeader();

    // Getters
    virtual double getObservable(int iObs);
    virtual std::string getObservableName() { return m_observableName; }
    virtual int getHeaderWidth() { return m_headerWidth; }

    // Setters
    virtual void reset();
    virtual void setLatticeSize(int latticeSize);
    virtual void setLatticeSpacing(double a) { m_a = a; }
    void setN(unsigned int *N);
    virtual void storeFlow(bool storeFlowObservable);// { m_storeFlowObservable = storeFlowObservable; } flytt denne til initialisering
};


#endif // CORRELATOR_H
