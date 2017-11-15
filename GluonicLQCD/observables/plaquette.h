#ifndef PLAQUETTE_H
#define PLAQUETTE_H

#include "correlator.h"

class Plaquette : public Correlator
{
private:
    double m_multiplicationFactor;
    SU3 P;
    SU3 PTemp;
    const std::string m_observableName = "Plaquette";
public:
    Plaquette(bool storeFlowObservable);
    ~Plaquette();
    void calculate(Links *lattice, int iObs);
    void calculate(SU3 *plaquetteStaples, int iObs);

    // Printers
    void printStatistics();
    // Setters
    void setLatticeSize(int latticeSize);
    void setPlaquetteStaples(SU3 *plaquetteStaples);
    virtual std::string getObservableName() { return m_observableName; }
};

#endif // PLAQUETTE_H
