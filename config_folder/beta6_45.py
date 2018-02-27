{   "bin_fn"                    : "build/GluonicLQCD",
    "runName"                   : "prodRunBeta6_45",
    "N"                         : 48,
    "NT"                        : 96,
    "subDims"                   : [12, 12, 12, 6],
    "beta"                      : 6.45,
    "NCf"                       : 250,
    "NCor"                      : 1600,
    "NTherm"                    : 20000,
    "NFlows"                    : 1000,
    "NUpdates"                  : 30,
    "storeCfgs"                 : True,
    "storeThermCfgs"            : False,
    "verboseRun"                : False,
    "hotStart"                  : False,
    "expFunc"                   : "morningstar", # options: luscher, taylor2, taylor4
    "observables"               : ["plaq"], # Optional: topologicalCharge, energyDensity
    "flowObservables"           : ["plaq","topc","energy","topct"], # Optional: topologicalCharge, energyDensity
    "uTest"                     : False,
    "uTestVerbose"              : False,
    "SU3Eps"                    : 0.24,
    "flowEpsilon"               : 0.02,
    "metropolisSeed"            : 0,
    "randomMatrixSeed"          : 0,
    "threads"                   : 512,
    "cpu_approx_runtime_hr"     : 360, # FLOW TIME: 159.4 hours +/- 12.830791 hours 
    "cpu_approx_runtime_min"    : 0,
    "cpu_memory"                : 3800,
    "account_name"              : "nn2977k"}
