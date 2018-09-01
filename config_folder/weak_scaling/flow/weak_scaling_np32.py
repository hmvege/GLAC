{
    "bin_fn": "build/GluonicLQCD", 
    "runName": "weak_scaling_np32_flow", 
    "N": 8, 
    "NT": 16, 
    "subDims": [
        4, 
        4, 
        4, 
        4
    ], 
    "beta": 6.0, 
    "NCf": 0, 
    "NCor": 0, 
    "NTherm": 0, 
    "NFlows": 1000, 
    "NUpdates": 0, 
    "storeCfgs": True, 
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
    "threads": 32, 
    "cpu_approx_runtime_hr": 3, 
    "cpu_approx_runtime_min": 0, 
    "cpu_memory": 3800, 
    "account_name": "nn2977k"
}