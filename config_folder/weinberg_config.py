# Run by:
#   python createJobs.py load weinberg_config.py -bf {base directory}
# 
# -bf {base_directory} is the full path leading up to where to output directory is situated.
# On my account, that would be /mnt/ls15/scratch/users/f0003959/ aka my scratch directory
# 
# This will submit a job to laconia with these settings.
#
# If you help, type -h, and you will get additional command line arguments.
#
# If you want to perform a dryrun(i.e. no permanent change is performed), run:
#   python createJobs.py --dryrun load weinberg_config.py -bf ${base directory}
# This will then print out every detail of the run, without submitting the job. Good if you 
# need to double check the settings.

{   
    "bin_fn"                    : "build/GluonicLQCD",
    "runName"                   : "weinbergConfigFlow1",
    "N"                         : 24,
    "NT"                        : 48,
    "subDims"                   : [6, 6, 6, 6],
    "beta"                      : 6.0,
    "NFlows"                    : 1000,
    "flowObservables"           : ["weinberg"],
    "flowEpsilon"               : 0.01,
    "inputFolder"               : "input", # Where to load the configurations from, can be absolute folder or relative folder if you wish.
    "outputFolder"              : "output", # Output directory. If not base-folder is given, this will be the same as the folder where you submit from
    "threads"                   : 512,
    "cpu_approx_runtime_hr"     : 65,
    "cpu_approx_runtime_min"    : 0,
    "cpu_memory"                : 3800,
    "account_name"              : "nn2977k",
}
