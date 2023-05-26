/*!
 * \class ObservableStorer
 *
 * \brief A class for storing the observable data.
 *
 * Contains methods for printing statistics, writing to file, setting observables and reseting ect.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef OBSERVABLESTORER_H
#define OBSERVABLESTORER_H

#include <string>
#include <vector>

class ObservableStorer
{
private:
    // Observable name
    std::string m_observableName;
    // Bool to store if we are to normalize the data by number of processors
    bool m_normalizeObservableByProcessor = false;

    // Observable data
    unsigned long m_NObs;
    double m_averagedObservable = 0;
    double m_varianceObservable = 0;
    double m_averagedObservableSquared = 0;
    double m_stdObservable = 0;
    std::vector<double> m_observables;
    std::vector<double> m_observablesSquared;
public:
    ObservableStorer(const unsigned long int NSize);
    ~ObservableStorer();

    // Accessor for the observable
    double &operator[](const unsigned long int iObs) { return m_observables.at(iObs); }
//    double &operator[](unsigned long int iObs) { return m_observables[iObs]; }

    // Runs statistics, perhaps create its own class? But that increases overhead, so maybe not
    void gatherResults();
    void runStatistics();

    // Printers
    void printStatistics();

    // File writers
    void writeObservableToFile(const double acceptanceRatio);
    void writeFlowObservableToFile(const unsigned long int configNumber);

    // Getters
    std::vector<double> getObservableArray() { return m_observables; }
    double getObservable(const unsigned long int iObs) { return m_observables[iObs]; }
    std::string getObservableName() { return m_observableName; }

    // Setters
    void setObservableName(const std::string &observableName) { m_observableName = observableName; }
    void setNormalizeObservableByProcessor(const bool norm) { m_normalizeObservableByProcessor = norm; }
    void reset();
};

#endif // OBSERVABLESTORER_H
