#include "configloader.h"
#include <fstream>
#include "lib/json.hpp"
#include "parameters.h"
#include "parallelization/parallel.h"


using nlohmann::json;

namespace ConfigLoader {
    namespace {
        // Namespace unaccessible outside of the ConfigLoader namespace.
        void setObservable(json iterable, bool flow) {
            std::vector<std::string> observableVector;
            for (unsigned int i = 0; i < iterable.size(); i++) {
                observableVector.push_back(iterable[i]);
            }
            if (flow) {
                Parameters::setFlowObservablesList(observableVector);
            } else {
                Parameters::setObservableList(observableVector);
            }
        }

        void setFieldConfigurations(json iterable) {
            std::vector<std::string> fieldConfigFileNames;
            for (unsigned int i = 0; i < iterable.size(); i++) {
                fieldConfigFileNames.push_back(iterable[i]);
            }
            Parameters::setFieldConfigurationFileNames(fieldConfigFileNames);
        }
    }

    void load(std::string jsonFileName)
    {
        std::ifstream i(jsonFileName.c_str());
        json j;
        i >> j;

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
        Parameters::setStoreConfigurations(bool(j["storeConfigurations"]));
        Parameters::setStoreThermalizationObservables(bool(j["storeThermalizationObservables"]));

        // Human readable output related variables
        Parameters::setVerbose(bool(j["verbose"]));

        // Setup related variables
        Parameters::setFilePath(j["pwd"]);
        Parameters::setBatchName(j["batchName"]);
        Parameters::setHotStart(bool(j["hotStart"]));
        Parameters::setRSTHotStart(bool(j["RSTHotStart"]));

        // Sub dimension setting
        if (!j["subDims"].empty()) {
            std::vector<unsigned int> tempN = {j["subDims"][0],j["subDims"][1],j["subDims"][2],j["subDims"][3]};
            Parameters::setN(tempN);
            Parameters::setSubLatticePreset(true);
            Parallel::Communicator::setN(tempN);
            Parallel::Index::setN(tempN);
        }

        // Exp.func setting (for flow)
        Parameters::setExpFuncName(j["expFunc"]);

        // Observables setting
        setObservable(j["observables"], false);

        // Flow observables setting
        setObservable(j["flowObservables"], true);

        // Loads field configurations if they are present
        if (j["load_field_configs"]) {
            setFieldConfigurations(j["field_configs"]);
            Parameters::setLoadFieldConfigurations(bool(j["load_field_configs"]));
            Parameters::setLoadChromaConfigurations(bool(j["chroma_config"]));
        }

        // Sets a field configuration to load, and then run the metropolis algorithm from the loaded configuration that is assumed to be thermalized
        std::string config_to_load = j["load_config_and_run"];
        if (config_to_load.length() != 0) { // If this is not empty, I will load and run a configuration
            std::vector<std::string> tempVec;
            tempVec.push_back(j["load_config_and_run"]);
            setFieldConfigurations(tempVec);
            Parameters::setNCf(j["NCf"]);
            Parameters::setLoadConfigAndRun(true);
        }

        if (!j["config_start_number"].empty()) {
            Parameters::setConfigStartNumber(j["config_start_number"]);
        }

        // Unit testing related variables
        Parameters::setUnitTesting(j["unitTesting"]);
        Parameters::setUnitTestingVerbose(j["unitTestingVerbose"]);

        // Performance testing related variables
        Parameters::setPerformanceTesting(j["performanceTesting"]);
        Parameters::setNExpTests(j["NExpTests"]);
        Parameters::setNRandTests(j["NRandTests"]);
        Parameters::setNDerivaitveTests(j["NDerivativeTests"]);
        Parameters::setTaylorPolDegree(j["TaylorPolDegree"]);

        // Checking if we have provideda gauge field to test
        std::string test_gauge_field = j["uTestFieldGaugeInvarince"];
        if (test_gauge_field.length() != 0) {
            Parameters::setCheckFieldGaugeInvariance(true);
            Parameters::setGaugeFieldToCheck(j["uTestFieldGaugeInvarince"]);
        }
        // Data generation related variables
        Parameters::setSU3Eps(j["SU3Eps"]);
        Parameters::setFlowEpsilon(j["flowEpsilon"]);
        if (j["metropolisSeed"] != 0) {
            Parameters::setMetropolisSeed(j["metropolisSeed"]);
        } else {
            Parameters::setMetropolisSeed(std::time(nullptr) + double(Parallel::Communicator::getNumProc()) + double(Parallel::Communicator::getProcessRank()));
        }
        if (j["randomMatrixSeed"] != 0) {
            Parameters::setRandomMatrixSeed(j["randomMatrixSeed"]);
        } else {
            Parameters::setRandomMatrixSeed(std::time(nullptr) + double(Parallel::Communicator::getProcessRank()));
        }
    }
}
