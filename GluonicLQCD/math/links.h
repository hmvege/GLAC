#ifndef LINKS_H
#define LINKS_H

#include "matrices/su3.h"

struct Links
{
    Links();
    ~Links();
    // Link variables being the matrices
    SU3 U[4];
};

#endif // LINKS_H
