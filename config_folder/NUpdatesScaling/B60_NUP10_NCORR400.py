{   "bin_fn"                    : "build/GluonicLQCD",
    "runName"                   : "B60_NUP10_NCORR400",
    "N"                         : 16,
    "NT"                        : 32,
    "subDims"                   : [4, 4, 4, 4],
    "beta"                      : 6.0,
    "NCf"                       : 200,
    "NCor"                      : 400,
    "NTherm"                    : 0,
    "NFlows"                    : 250,
    "NUpdates"                  : 10,
    "storeCfgs"                 : False,
    "storeThermCfgs"            : False,
    "verboseRun"                : False,
    "hotStart"                  : False,
    "expFunc"                   : "morningstar", # options: luscher, taylor2, taylor4
    "observables"               : ["plaq"], # Optional: topologicalCharge, energyDensity
    "flowObservables"           : ["plaq","topc","energy"], # Optional: topologicalCharge, energyDensity
    "uTest"                     : False,
    "uTestVerbose"              : False,
    "SU3Eps"                    : 0.24,
    "flowEpsilon"               : 0.01,
    "metropolisSeed"            : 0,
    "randomMatrixSeed"          : 0,
    "threads"                   : 512,
    "cpu_approx_runtime_hr"     : 4,
    "cpu_approx_runtime_min"    : 0,
    "cpu_memory"                : 3800,
    "account_name"              : "nn2977k"}
