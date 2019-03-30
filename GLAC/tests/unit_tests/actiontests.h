#ifndef ACTIONTESTS_H
#define ACTIONTESTS_H

#include "testcore.h"

class ActionTests : public TestCore
{
private:
    // Tests for action
    bool testWilsonAction();
    bool testExplicitExpAction();
public:
    ActionTests();

    bool runActionTests();
};

#endif // ACTIONTESTS_H
