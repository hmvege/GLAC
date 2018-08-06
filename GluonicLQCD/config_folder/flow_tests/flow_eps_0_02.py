{   "bin_fn"                    : "build/GluonicLQCD",
    "runName"                   : "flow_eps_0_02", # 0.001 0.005 0.007 0.01 0.02 0.03 0.05 0.1 0.5
    "N"                         : 24,
    "NT"                        : 48,
    "subDims"                   : [6, 6, 6, 6],
    "beta"                      : 6.0,
    "NCf"                       : 1, # LOADS A CONFIG FROM feps test config generator!!
    "NCor"                      : 600,
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
    "flowEpsilon"               : 0.02, # 0.001 0.005 0.007 0.01 0.02 0.03 0.05 0.1 0.5
    "metropolisSeed"            : 0,
    "randomMatrixSeed"          : 0,
    "threads"                   : 512,
    "cpu_approx_runtime_hr"     : 2,
    "cpu_approx_runtime_min"    : 30,
    "cpu_memory"                : 3800,
    "account_name"              : "nn2977k"}
