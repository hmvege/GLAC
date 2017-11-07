#include "topologicalcharge.h"
#include "functions.h"
#include "clover.h"
#include <cmath>

TopologicalCharge::TopologicalCharge() : Correlator()
{
    m_multiplicationFactor = 1.0/(32*M_PI*M_PI);
    populateLC(); // Fills the levi civita vector
}

TopologicalCharge::~TopologicalCharge()
{

}

double TopologicalCharge::calculate(Links *lattice)
{
    /*
     * Function to be used when no clover is provided. SHOULD BE TESTED
     */
    Clover Clov;
    Clov.initializeIndexHandler(m_Index);
    Clov.setN(m_N);
    Clov.setLatticeSize(m_latticeSize);
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
                    for (unsigned int i = 0; i < m_leviCivita.size(); i++)
                    {
                        G1 = Clov.m_clovers[m_leviCivita[i].ci[0]];
                        G2 = Clov.m_clovers[m_leviCivita[i].ci[1]];
                        topCharge += traceSparseImagMultiplication(G1,G2)*m_leviCivita[i].sgn;
                    }
                }
            }
        }
    }
    return topCharge*m_multiplicationFactor;
}

double TopologicalCharge::calculate(SU3 *clovers)
{
    topCharge = 0;
    for (unsigned int i = 0; i < m_leviCivita.size(); i++)
    {
        G1 = clovers[m_leviCivita[i].ci[0]];
        G2 = clovers[m_leviCivita[i].ci[1]];
        topCharge += traceSparseImagMultiplication(G1,G2)*m_leviCivita[i].sgn;
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
