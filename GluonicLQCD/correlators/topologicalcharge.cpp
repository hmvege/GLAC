#include "topologicalcharge.h"
#include "functions.h"
#include <cmath>


TopologicalCharge::TopologicalCharge() : Correlator()
{
    m_multiplicationFactor = 1.0/(32*std::atan(1)*4*std::atan(1)*4);
    populateLC(); // Fills the levi civita vector
}

TopologicalCharge::~TopologicalCharge()
{

}

void TopologicalCharge::setClover(SU3 *clover)
{
    for (int i = 0; i < 12; i++) {
        m_clover[i] = clover[i];
    }
}

double TopologicalCharge::calculate()
{
//    for (int mu = 0; mu < 4; mu++) {
//        for (int nu = 0; nu < 4; nu++) {
//            for (int rho = 0; rho < 4; rho++) {
//                for (int sigma = 0; sigma < 4; sigma++) {

//                }
//            }
//        }
//    }
    double topCharge = 0;
    SU3 G1,G2;
    for (unsigned int i = 0; i < m_leviCivita.size(); i++)
    {
        G1 = m_clover[cloverIndex(m_leviCivita[i].lc[0],m_leviCivita[i].lc[1])];
        G2 = m_clover[cloverIndex(m_leviCivita[i].lc[2],m_leviCivita[i].lc[3])];
        topCharge += traceImagMultiplication(G1,G2)*m_leviCivita[i].sgn;
    }
    return topCharge*m_multiplicationFactor;
}

void TopologicalCharge::populateLC()
{
    int muNuOverCounter = 0;
    int rhoSigmaOverCounter = 0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            if (nu==mu) {
                muNuOverCounter++;
                continue;
            }
            for (int rho = 0; rho < 4; rho++) {
                if (rho==mu || rho==nu) continue;
                for (int sigma = 0; sigma < 4; sigma++) {
                    if (sigma==rho) {
                        rhoSigmaOverCounter++;
                        continue;
                        if (sigma==mu || sigma==nu) continue;
                    }
                    LeviCivita LC;
                    LC.setLC(mu, nu - muNuOverCounter, rho, sigma - rhoSigmaOverCounter);
                    LC.sgn = getLCSign(LC);
                    m_leviCivita.push_back(LC);
                }
            }
        }
    }
}

int TopologicalCharge::getLCSign(LeviCivita LC)
{
    int sign = 1;
    int x = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            x = (LC.lc[i] - LC.lc[j]);
            sign *= (x > 0) - (x < 0);
        }
    }
    return sign;
}
