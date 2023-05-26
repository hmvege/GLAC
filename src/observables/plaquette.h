/*!
 * \class Plaquette
 *
 * \brief Observable for calculating the plaquette.
 *
 * Plaqeutte definition:
 * \f[
 *      P = \frac{1}{16|\Lambda|} \sum_{n\in\Lambda} \sum_{\mu < \nu} \Re \mathrm{tr} P_{\mu\nu},
 * \f]
 * where
 * \f[
 *      P_{\mu\nu}=U_\mu(n) U_{\nu}(n+\hat{\mu}) U_{\mu}(n+\hat{\nu})^\dagger U_{\nu} (n)^\dagger
 * \f]
 * is the plquette.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
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
    Plaquette(const bool storeFlowObservable);
    ~Plaquette();
    void calculate(Lattice<SU3> *lattice, const unsigned int iObs);
    // Statistics getter
    void runStatistics();
    // Setters
    void setLatticeSize(const unsigned long int latticeSize);
    // Getters
    std::string getObservableName() { return m_observableName; }
    // Printers
    void printObservable(const unsigned int iObs);
    void printHeader();
};

#endif // PLAQUETTE_H
