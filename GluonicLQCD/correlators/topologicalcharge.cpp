#include "topologicalcharge.h"
#include "functions.h"
#include "clover.h"
#include <cmath>

TopologicalCharge::TopologicalCharge(double a) : Correlator()
{
    m_a = a;
    m_multiplicationFactor = m_a*m_a*m_a*m_a/(32*std::atan(1)*4*std::atan(1)*4);
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

double TopologicalCharge::calculate(Links *lattice)
{
    /*
     * Function to be used when no clover is provided. SHOULD BE TESTED
     */
    Clover Clov;
    topCharge = 0;
    for (unsigned int i = 0; i < m_N[0]; i++) { // x
        for (unsigned int j = 0; j < m_N[1]; j++) { // y
            for (unsigned int k = 0; k < m_N[2]; k++) { // z
                for (unsigned int l = 0; l < m_N[3]; l++) { // t
                    m_position[0] = i;
                    m_position[1] = j;
                    m_position[2] = k;
                    m_position[3] = l;
                    Clov.calculateClover(lattice,i,j,k,l);
                    setClover(Clov.m_clovers);
                    for (unsigned int i = 0; i < m_leviCivita.size(); i++)
                    {
                        G1 = m_clover[cloverIndex(m_leviCivita[i].lc[0],m_leviCivita[i].lc[1])];
                        G2 = m_clover[cloverIndex(m_leviCivita[i].lc[2],m_leviCivita[i].lc[3])];
                        topCharge += traceImagMultiplication(G1,G2)*m_leviCivita[i].sgn;
                    }
                }
            }
        }
    }
    return topCharge*m_multiplicationFactor*m_a*m_a*m_a*m_a;
}

double TopologicalCharge::calculate()
{
    topCharge = 0;
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
                muNuOverCounter++; // Acounts for overcounting
                continue;
            }

            for (int rho = 0; rho < 4; rho++) {
                if (rho==mu || rho==nu) {
                    rhoSigmaOverCounter++; // Acounts for overcounting
                    continue;
                }

                for (int sigma = 0; sigma < 4; sigma++) {

                    if (sigma==rho) {
                        rhoSigmaOverCounter++; // Acounts for overcounting
                        continue;
                    }
                    if (sigma==mu || sigma==nu) continue;
                    LeviCivita LC;
                    LC.setLC(mu, nu, rho, sigma);
                    LC.sgn = getLCSign(LC);
                    LC.setCI(cloverIndex(mu,nu - muNuOverCounter),cloverIndex(rho,sigma - rhoSigmaOverCounter));
                    m_leviCivita.push_back(LC);
                }
            }
            rhoSigmaOverCounter = 0;
        }
    }
    if (m_leviCivita.size() != 24)
    {
        cout << "Error: number of levi civita combinations is " << m_leviCivita.size() << ", not 24!" << endl;
        for (unsigned int i = 0; i < m_leviCivita.size(); i++) {
            cout << m_leviCivita[i];
            cout << "   Clover index: " << cloverIndex(m_leviCivita[i].lc[0],m_leviCivita[i].lc[1]) << " " << cloverIndex(m_leviCivita[i].lc[2],m_leviCivita[i].lc[3]) << endl;
        }
        exit(1);
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
