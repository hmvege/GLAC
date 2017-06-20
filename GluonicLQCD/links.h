#ifndef LINKS_H
#define LINKS_H

#include "su3.h"

class Links
{
public:
    Links();
    ~Links();
    // Link variables being the matrices
    SU3 * U;
};

#endif // LINKS_H
