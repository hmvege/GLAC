#ifndef TOPOLOGICALCHARGE_H
#define TOPOLOGICALCHARGE_H

#include "correlator.h"

class TopologicalCharge : public Correlator
{
private:
    struct LeviCivita {
        /*
         * Small struct for containing the levi civita symbols.
         */
        int lc[4];
        int sgn;
        void setLC(int mu, int nu, int rho, int sigma) {
            lc[0] = mu;
            lc[1] = nu;
            lc[2] = rho;
            lc[3] = sigma;
        }
//        int getMu() { return lc[0]; }
//        int getNu() { return lc[1]; }
//        int getRho() { return lc[2]; }
//        int getSigma() { return lc[3]; }
        friend std::ostream& operator<<(std::ostream& os, const LeviCivita& a) {
            os << "("<< a.lc[0] << " " << a.lc[1] << " " << a.lc[2] << " " << a.lc[3] << ") sign: " << a.sgn;
            return os;
        }
    };

    SU3 m_clover[12];
    double m_multiplicationFactor;
    std::vector<LeviCivita> m_leviCivita;
//    LeviCivita m_leviCevita[24];
//    LeviCivita m_leviCevitaOdd[12];

    void populateLC();
    int getLCSign(LeviCivita LC);
public:
    TopologicalCharge();
    ~TopologicalCharge();

    double calculate();
    void setClover(SU3 *clover);
};

#endif // TOPOLOGICALCHARGE_H
