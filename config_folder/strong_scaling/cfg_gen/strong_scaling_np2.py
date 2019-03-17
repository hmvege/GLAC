{
    "bin_fn": "build/GluonicLQCD", 
    "runName": "strong_scaling_np2_gen", 
    "N": 32, 
    "NT": 16, 
    "subDims": [
        32, 
        32, 
        16, 
        16
    ], 
    "beta": 6.0, 
    "NCf": 1, 
    "NCor": 600, 
    "NTherm": 0, 
    "NFlows": 0, 
    "NUpdates": 30, 
    "storeCfgs": False, 
    "storeThermCfgs": False, 
    "verboseRun": False, 
    "hotStart": False, 
    "expFunc": "morningstar", 
    "observables": [
        "plaq"
    ], 
    "flowObservables": [
        "plaq", 
        "topc", 
        "energy", 
        "topct"
    ], 
    "uTest": False, 
    "uTestVerbose": False, 
    "SU3Eps": 0.24, 
    "flowEpsilon": 0.01, 
    "metropolisSeed": 0, 
    "randomMatrixSeed": 0, 
    "threads": 2, 
    "cpu_approx_runtime_hr": 64, 
    "cpu_approx_runtime_min": 0, 
    "cpu_memory": 3800, 
    "account_name": "nn2977k"
}