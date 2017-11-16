{   "bin_fn"                    : "build/GluonicLQCD",
    "runName"                   : "sizeScaling32cube64",
    "N"                         : 32,
    "NT"                        : 64,
    "subDims"                   : [16, 16, 16, 8],
    "beta"                      : 6.0,
    "NCf"                       : 1,
    "NCor"                      : 2000,
    "NTherm"                    : 0,
    "NFlows"                    : 0,
    "NUpdates"                  : 10,
    "storeCfgs"                 : 0,
    "storeThermCfgs"            : 0,
    "verboseRun"                : 1,
    "hotStart"                  : 0,
    "expFunc"                   : "morningstar", # options: luscher, taylor2, taylor4
    "observables"               : ["plaquette"], # Optional: topologicalCharge, energyDensity
    "flowObservables"           : ["plaquette"], # Optional: topologicalCharge, energyDensity
    "uTest"                     : 0,
    "uTestVerbose"              : 0,
    "SU3Eps"                    : 0.24,
    "flowEpsilon"               : 0.01,
    "metropolisSeed"            : 0,
    "randomMatrixSeed"          : 0,
    "threads"                   : 64,
    "cpu_approx_runtime_hr"     : 2,
    "cpu_approx_runtime_min"    : 0,
    "cpu_memory"                : 3800,
    "account_name"              : "nn2977k"}