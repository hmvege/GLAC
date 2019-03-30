#include "actiontests.h"
#include "actions/actions.h"
#include "parallelization/communicator.h"

ActionTests::ActionTests()
{

}

bool ActionTests::testWilsonAction()
{
    bool passed = true;

    WilsonGaugeAction S;

    Lattice<SU3> * L = new Lattice<SU3>[4];
    for (int i = 0; i < 4; i++) L[i].allocate(m_dim);

    for (int mu = 0; mu < 4; mu++)
    {
        for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
        {
            L[mu][iSite].identity();
        }
    }

    // Tests the delta action, as used in cfg generation. The staple is indirectly tested.
    SU3 updateLink;
    updateLink.identity();
    updateLink *= 2;

    double dS = 0;

    for (unsigned int x = 0; x < L->m_dim[0]; x++) {
        for (unsigned int y = 0; y < L->m_dim[1]; y++) {
            for (unsigned int z = 0; z < L->m_dim[2]; z++) {
                for (unsigned int t = 0; t < L->m_dim[3]; t++) {
                    for (int mu = 0; mu < 4; mu++) {
                        S.computeStaple(L, x, y, z, t, mu);
                        dS = S.getDeltaAction(L[mu][Parallel::Index::getIndex(x,y,z,t)], updateLink);
                        if (fabs(dS + 36) > 1e-15) {
                            cout << "Error in action: " << dS << endl;
                            passed = false;
                        }
                    }
                }
            }
        }
    }

//    // Tests the action derivative.
//    L[0][0].print();
//    Lattice<SU3> L1 = S.getActionDerivative(L, 0);
//    L1[0].print();

    delete [] L;

    return passed;
}

bool ActionTests::testExplicitExpAction()
{
    bool passed = true;

    return passed;
}

// Action tests
bool ActionTests::runActionTests()
{
    if (m_processRank == 0 && m_verbose) {
        printf("Running gauge action tests.\n");
    }

    bool passed = (testWilsonAction());// & testExplicitExpAction());

    if (passed) {
        if (m_processRank == 0) cout << "PASSED: gauge action." << endl;
    } else {
        if (m_processRank == 0) cout << "FAILED: gauge action." << endl;
    }
    Parallel::Communicator::setBarrier();

    return passed;
}
