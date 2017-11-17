#include "configloader.h"

using nlohmann::json;

void ConfigLoader::load(std::string jsonFileName)
{
    std::ifstream i(jsonFileName.c_str());
    json j;
    i >> j;
    if (Parallel::Communicator::getProcessRank() == 0) std::cout << std::setw(4) << j << std::endl;

    // Lattice related run variables
    Parameters::setNSpatial(j["NSpatial"]);
    Parameters::setNTemporal(j["NTemporal"]);
    Parallel::Index::setNTot(j["NSpatial"],j["NTemporal"]);
    Parameters::setBeta(j["beta"]);
    Parameters::setNCf(j["NCf"]);
    Parameters::setNCor(j["NCor"]);
    Parameters::setNTherm(j["NTherm"]);
    Parameters::setNFlows(j["NFlows"]);
    Parameters::setNUpdates(j["NUpdates"]);
    // Data storage related variables
    Parameters::setOutputFolder(j["outputFolder"]);
    Parameters::setInputFolder(j["inputFolder"]);
    Parameters::setStoreConfigurations(j["storeConfigurations"]);
    Parameters::setStoreThermalizationObservables(j["storeThermalizationObservables"]);
    // Human readable output related variables
    Parameters::setVerbose(j["verbose"]);
    // Setup related variables
    Parameters::setFilePath(j["pwd"]);
    Parameters::setBatchName(j["batchName"]);
    Parameters::setHotStart(j["hotStart"]);
    // Sub dimension setting
    if (!j["subDims"].empty()) {
        unsigned int tempN[4] = {j["subDims"][0],j["subDims"][1],j["subDims"][2],j["subDims"][3]};
        Parameters::setN(tempN);
        Parallel::Index::setN(tempN);
        Parameters::setSubLatticePreset(true);
    }
    // Exp.func setting (for flow)

    // Observables setting

    // Flow observables setting

    // Testing related variables
    Parameters::setUnitTesting(j["unitTesting"]);
    Parameters::setUnitTestingVerbose(j["unitTestingVerbose"]);
    // Data generation related variables
    Parameters::setSU3Eps(j["SU3Eps"]);
    Parameters::setFlowEpsilon(j["flowEpsilon"]);
    Parameters::setMetropolisSeed(j["metropolisSeed"]);
    Parameters::setRandomMatrixSeed(j["randomMatrixSeed"]);
}
