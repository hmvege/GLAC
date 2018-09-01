{
    "bin_fn": "build/GluonicLQCD", 
    "runName": "weak_scaling_np2_io", 
    "N": 4, 
    "NT": 8, 
    "subDims": [
        4, 
        4, 
        4, 
        4
    ], 
    "beta": 6.0, 
    "NCf": 10, 
    "NCor": 1, 
    "NTherm": 0, 
    "NFlows": 0, 
    "NUpdates": 1, 
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
    "threads": 2, 
    "cpu_approx_runtime_hr": 2, 
    "cpu_approx_runtime_min": 0, 
    "cpu_memory": 3800, 
    "account_name": "nn2977k"
}