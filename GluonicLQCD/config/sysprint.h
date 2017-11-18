#ifndef SYSPRINT_H
#define SYSPRINT_H

#include <iomanip>
#include <iostream>
#include "parameters.h"
#include "parallelization/communicator.h"

class SysPrint
{
private:
    static std::string getTrueOrFalseString(bool i);
    static std::string getListString(std::vector<std::string> observableList);
public:
    SysPrint();

    static void printSystemInfo();
    static void printLine();
};

#endif // SYSPRINT_H
