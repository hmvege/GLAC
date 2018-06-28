{   "bin_fn"                    : "build/GluonicLQCD", # INTENDED FOR SMAUG
    "runName"                   : "topc_modes_run_beta60",
    "N"                         : 16,
    "NT"                        : 32,
    "subDims"                   : [8, 8, 8, 4],
    "beta"                      : 6.0,
    "NCf"                       : 2500,
    "NCor"                      : 200,
    "NTherm"                    : 10000,
    "NFlows"                    : 1000,
    "NUpdates"                  : 10,
    "storeCfgs"                 : True,
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
    "threads"                   : 64,
    "cpu_approx_runtime_hr"     : 200, # FLOW TIME: 70.0 hours
    "cpu_approx_runtime_min"    : 0,
    "cpu_memory"                : 2000,
    "account_name"              : "nn2977k"}
