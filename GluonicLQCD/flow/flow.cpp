#include "flow.h"
#include "links.h"
#include "parallelization/indexorganiser.h"

/*
 * A class for performing Wilson flow on a gauge field configuration.
 */

Flow::flow(int *N)
{
    for (int i = 0; i < 4; i++) m_N[i] = N[i];
}

void Flow::flowGaugeField(int NFlows, Links *lattice)
{

    for (int i = 0; i < NFlows; i++) {
        runFlow(lattice);
    }
}

void Flow::runFlow(Links *lattice)
{
    for (unsigned int x = 0; x < m_N[0]; x++) {
        for (unsigned int y = 0; y < m_N[1]; y++) {
            for (unsigned int z = 0; z < m_N[2]; z++) {
                for (unsigned int t = 0; t < m_N[3]; t++) {
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        smearLink();
                    }
                }
            }
        }
    }
}

void Flow::smearLink()
{

}

void Flow::setIndexHandler(IndexOrganiser *Index)
{
    m_Index = Index;
}
