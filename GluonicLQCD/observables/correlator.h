#ifndef CORRELATOR_H
#define CORRELATOR_H

#include "math/lattice.h"
#include "tools/observablestorer.h"
#include <vector>

class Correlator
{
protected:
    // Lattice dimension array
    std::vector<unsigned int> m_N;

    // (sub)lattice size
    double m_latticeSize;

    // Lattice spacing
    double m_a;

    // Observable name
    const int m_headerWidth = 20;
    const std::string m_observableName = "Correlator";
    const std::string m_observableNameCompact = "correlator";

    // Creates a object that store the observable
    ObservableStorer * m_observable = nullptr;

    // Object that holds info of whether or not we are storing flow observables
    bool m_storeFlowObservable = false;
public:
    Correlator(bool storeFlowObservable);
    Correlator();
    virtual ~Correlator();
    virtual void calculate(Lattice<SU3> *lattice, int iObs);
    virtual void writeObservableToFile(double acceptanceRatio);
    virtual void writeFlowObservablesToFile(int iFlow);
    virtual void runStatistics();

    // Printers
    virtual void printObservable(int iObs);
    virtual void printHeader();
    virtual void printStatistics();

    // Observable copyer
    virtual void copyObservable(int iObs, std::vector<double> obs);
    virtual std::vector<double> getObservablesVector(int iObs);

    // Getters
    virtual double getObservable(int iObs);
    virtual std::string getObservableName() { return m_observableName; }
    virtual int getHeaderWidth() { return m_headerWidth; }

    // Setters
    virtual void reset();
    virtual void setObservable(int iObs, double obs);
    virtual void storeFlow(bool storeFlowObservable);
};


#endif // CORRELATOR_H
