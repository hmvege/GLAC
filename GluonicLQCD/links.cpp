#include "links.h"

/*
 * For storing a link variable
 */

Links::Links()
{
    U = new SU3[4];
}

Links::~Links()
{
    delete [] U;
}
