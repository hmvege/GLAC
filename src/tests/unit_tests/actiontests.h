#ifndef ACTIONTESTS_H
#define ACTIONTESTS_H

#include "testcore.h"
#include "actions/action.h"

class ActionTests : public TestCore
{
private:
    // Tests for action
    bool testWilsonAction();
    bool testExplicitExpAction();

    bool testActionDerivative(Action &S);
    bool testActionDeltaS(Action &S);
public:
    ActionTests();

    bool runActionTests();
};

#endif // ACTIONTESTS_H
