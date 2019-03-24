/*!
 * \namespace ConfigLoader
 *
 * \brief json configuration loader.
 *
 * Based upon json.hpp. Loads a json configurations, and sets the parameters in the Parameters class.
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef CONFIGLOADER_H
#define CONFIGLOADER_H

#include <string>

namespace ConfigLoader
{
    void load(std::string jsonFileName);
}

#endif // CONFIGLOADER_H
