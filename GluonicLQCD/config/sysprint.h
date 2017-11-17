#ifndef SYSPRINT_H
#define SYSPRINT_H

#include <iomanip>
#include <iostream>
#include "parameters.h"
#include "parallelization/communicator.h"

class SysPrint
{
public:
    SysPrint();

    static void printSystemInfo();
    inline static void printLine();
};

#endif // SYSPRINT_H
