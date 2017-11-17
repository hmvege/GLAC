#ifndef CONFIGLOADER_H
#define CONFIGLOADER_H

#include <fstream>
#include "lib/json.hpp"
#include "parameters.h"
#include "parallelization/communicator.h"
#include "parallelization/index.h"

namespace ConfigLoader
{
    void load(std::string jsonFileName);
    namespace {
        // Namespace unaccessible outside of the ConfigLoader namespace.
    }
}

#endif // CONFIGLOADER_H
