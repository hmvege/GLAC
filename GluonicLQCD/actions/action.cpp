#include "action.h"
#include "links.h"
#include "matrices/su3.h"
#include "functions.h"

// TEMP
#include <iostream>
using std::cout;
using std::endl;

Action::Action(int N, int N_T)
{
    m_N = N;
    m_N_T = N_T;
}

Action::~Action()
{

}

double Action::getDeltaAction(Links * lattice, SU3 U, int i, int j, int k, int l, int mu)
{
    cout << "In Action::getDeltaAction: If you are seeing this, something is wrong!" << endl;
    exit(1);
    return 1.0;
}

void Action::computeStaple(Links *lattice, int i, int j, int k, int l, int mu)
{
    cout << "In Action::computeStaple: If you are seeing this, something is wrong!" << endl;
    exit(1);
}
