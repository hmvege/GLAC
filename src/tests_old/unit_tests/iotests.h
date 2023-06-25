#ifndef IOTESTS_H
#define IOTESTS_H

#include "testcore.h"

class IOTests : public TestCore
{
private:
    // Test for IO
    bool testIOLatticeWriteRead();
    bool testIOWriteDoubles();
public:
    IOTests();

    bool runIOTests();
};

#endif // IOTESTS_H
