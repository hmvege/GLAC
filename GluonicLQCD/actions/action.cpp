#include "action.h"

Action::Action()
{
    m_N.resize(4);
    for (int i = 0; i < 4; i++) {
        m_N[i] = 0;
    }
    m_position = std::vector<int>(4,0);
}

Action::~Action()
{
}

double Action::getDeltaAction(SU3 U, SU3 UPrime)
{
    cout << "In Action::getDeltaAction: If you are seeing this, something is wrong!" << endl;
    exit(1);
    U = UPrime;
    return 1.0;
}

void Action::computeStaple(Lattice<SU3> *lattice, unsigned int i, unsigned int j, unsigned int k, unsigned int l, int mu)
{
    cout << "In Action::computeStaple: If you are seeing this, something is wrong!" << endl;
    exit(1);
    lattice[mu][i+j+k+l].print();
}

void Action::setN(std::vector<unsigned int> N)
{
    m_N = N;
}

Lattice<SU3> Action::getActionDerivative(Lattice<SU3> *lattice, int mu)
{
    cout << "In Action::getActionDerivative: If you are seeing this, something is wrong!" << endl;
    exit(1);
    lattice[mu][0].print();
}
