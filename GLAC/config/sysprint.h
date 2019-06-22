/*!
 * \class SysPrint
 *
 * \brief Class for printing run information.
 *
 * Holds all of the parameters passed on from the json file. Stores parameters as static, thus being accessible from everywhere.
 *
 * \todo Make this class into a class that prints instead of having to do if rank=0 ect. Maybe move to IO as well.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef SYSPRINT_H
#define SYSPRINT_H

#include "parameters.h"
#include "parallelization/parallel.h"

class SysPrint
{
private:
    static std::string getTrueOrFalseString(const bool i);
    static std::string getListString(const std::vector<std::string> &observableList);
public:
    SysPrint();

    static void printSystemInfo();
    static void printLine();
};

#endif // SYSPRINT_H
