#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"

class Plaquette : public Correlator
{
private:
    double m_multiplicationFactor;
    double m_tempObservable = 0;
    // Initializes temporary sample storer
    Lattice<SU3> m_temp;
    // Observable names, human readability and for io
    const std::string m_observableName = "Plaquette";
    const std::string m_observableNameCompact = "plaq";
public:
    Plaquette(bool storeFlowObservable);
    ~Plaquette();
    void calculate(Lattice<SU3> *lattice, unsigned int iObs);
    // Statistics getter
    void runStatistics();
    // Setters
    void setLatticeSize(unsigned long int latticeSize);
    // Getters
    std::string getObservableName() { return m_observableName; }
    // Printers
    void printObservable(unsigned int iObs);
    void printHeader();
};

#endif // PLAQUETTE_H
