{   "bin_fn"                    : "build/GluonicLQCD",
    "runName"                   : "distribution_test_eps020",
    "N"                         : 6,
    "NT"                        : 12,
    "subDims"                   : [3, 3, 3, 3],
    "beta"                      : 6.0,
    "NCf"                       : 2000,
    "NCor"                      : 200,
    "NTherm"                    : 1000,
    "NFlows"                    : 1000,
    "NUpdates"                  : 20,
    "storeCfgs"                 : False,
    "storeThermCfgs"            : True,
    "verboseRun"                : False,
    "hotStart"                  : False,
    "expFunc"                   : "morningstar", # options: luscher, taylor2, taylor4
    "observables"               : ["plaq"], # Optional: topologicalCharge, energyDensity
    "flowObservables"           : ["plaq", "topc", "energy", "topct"], # Optional: topologicalCharge, energyDensity
    "uTest"                     : False,
    "uTestVerbose"              : False,
    "SU3Eps"                    : 0.20,
    "flowEpsilon"               : 0.01,
    "metropolisSeed"            : 0,
    "randomMatrixSeed"          : 0,
    "threads"                   : 32,
    "cpu_approx_runtime_hr"     : 100, # FLOW TIME: 38.3 hours
    "cpu_approx_runtime_min"    : 0,
    "cpu_memory"                : 3800,
    "account_name"              : "nn2977k"}
