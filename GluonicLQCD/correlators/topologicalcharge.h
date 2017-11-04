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
        friend std::ostream& operator<<(std::ostream& os, const LeviCivita& a) {
            os << "("<< a.lc[0] << " " << a.lc[1] << " " << a.lc[2] << " " << a.lc[3] << ") sign: " << a.sgn;
            return os;
        }
    };

    double topCharge;
    double m_a; // Lattice spacing
    SU3 G1,G2;
    SU3 m_clover[12];
    double m_multiplicationFactor;
    std::vector<LeviCivita> m_leviCivita;

    void populateLC();
    int getLCSign(LeviCivita LC);
public:
    TopologicalCharge(double a);
    ~TopologicalCharge();

    double calculate();
    double calculate(Links *lattice);
    void setClover(SU3 *clover);
};

#endif // TOPOLOGICALCHARGE_H
