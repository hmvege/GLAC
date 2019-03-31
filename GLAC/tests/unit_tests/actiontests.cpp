#include "actiontests.h"
#include "actions/actions.h"
#include "parallelization/communicator.h"

ActionTests::ActionTests()
{

}

bool ActionTests::testActionDerivative(Action &S)
{
    bool passed = true;

    Lattice<SU3> * L = new Lattice<SU3>[4];
    Lattice<SU3> * LNew = new Lattice<SU3>[4];
    for (int i = 0; i < 4; i++) L[i].allocate(m_dim);
    for (int i = 0; i < 4; i++) LNew[i].allocate(m_dim);


    for (int mu = 0; mu < 4; mu++)
    {
        L[mu].identity();
        LNew[mu].zeros();
    }

    // Tests the action derivative as used in flow.
    for (int mu = 0; mu < 4; mu++) {
        LNew[mu] = S.getActionDerivative(L, mu);
    }


//    SU3 updateLink;
//    updateLink.identity();
//    updateLink *= 2;

//    double dS = 0;

//    for (unsigned int x = 0; x < L->m_dim[0]; x++) {
//        for (unsigned int y = 0; y < L->m_dim[1]; y++) {
//            for (unsigned int z = 0; z < L->m_dim[2]; z++) {
//                for (unsigned int t = 0; t < L->m_dim[3]; t++) {
//                    for (int mu = 0; mu < 4; mu++) {
//                        S.computeStaple(L, x, y, z, t, mu);
//                        dS = S.getDeltaAction(L[mu][Parallel::Index::getIndex(x,y,z,t)], updateLink);

//                        if (fabs(dS + 36) > 1e-15) {
//                            cout << "Error in action(should equal 36): " << dS << endl;
//                            passed = false;
//                            x = L->m_dim[0];
//                            y = L->m_dim[1];
//                            z = L->m_dim[2];
//                            t = L->m_dim[3];
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//    }

    delete [] L;
    delete [] LNew;

    return passed;
}

bool ActionTests::testActionDeltaS(Action &S)
{
//    double getDeltaAction(SU3 U, SU3 UPrime);
//    void computeStaple(Lattice<SU3> *lattice, int i, int j, int k, int l, int mu);
    bool passed = true;

    Lattice<SU3> * L = new Lattice<SU3>[4];
    for (int i = 0; i < 4; i++) L[i].allocate(m_dim);

    for (int mu = 0; mu < 4; mu++)
    {
        L[mu].identity();
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
                            cout << "Error in action deltaS (should equal 36): " << dS << endl;
                            passed = false;
                            x = L->m_dim[0];
                            y = L->m_dim[1];
                            z = L->m_dim[2];
                            t = L->m_dim[3];
                            break;
                        }
                    }
                }
            }
        }
    }

    delete [] L;

    return passed;

}

bool ActionTests::testWilsonAction()
{
    bool passed;

    WilsonGaugeAction S;

    if (m_verbose && m_processRank == 0) {
        cout << "Testing WilsonGaugeAction." << endl;
    }
    passed = testActionDeltaS(S);
//    passed = testActionDerivative(S);

    return passed;
}

bool ActionTests::testExplicitExpAction()
{
    bool passed = true;

    WilsonExplicitExp S;

    if (m_verbose && m_processRank == 0) {
        cout << "Testing WilsonExplicitExp." << endl;
    }
    passed = testActionDeltaS(S);
//    passed = testActionDerivative(S);

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
