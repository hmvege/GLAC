#ifndef OBSERVABLETESTS_H
#define OBSERVABLETESTS_H

#include "testcore.h"

class ObservableTests : public TestCore
{
private:
    // Tests for observables
    bool testTopCharge();
public:
    ObservableTests();

    bool runObservableTests();
};

#endif // OBSERVABLETESTS_H
