#!/usr/bin/env python2

import os
import subprocess
import time
import sys
import argparse
import json
import ast
import shutil
import re

# TODO: implement a way of checking that the .bin file exists in case of lcfgr

AVAILABLE_OBSERVABLES = ["plaq", "topc", "energy", "topct", "weinberg", "weinbergt"]
# AVAILABLE_SCALAR_FIELDS = ["plaq", "topc", "energy", "weinberg"]
AVAILABLE_SCALAR_FIELDS = ["energyTopcFieldDensity"]
AVAILABLE_EXP_FUNCS = ["morningstar", "luscher", "taylor2", "taylor4"]
AVAILABLE_ACTIONS = ["luscher", "wilson"]
AVAILABLE_HPC_SYSTEMS = ["slurm", "torque", "local"]

def get_arg_max_index(N):
    """For getting the maximum index of an list."""
    val = N[0]
    index = 0
    for i in xrange(4):
        if N[i] > val:
            val = N[i]
            index = i
    return index

def create_square(numprocs, NSpatial, NTemporal):
    """
    Create a square sub lattice.

    Args:
        numprocs: integer, number of processors.
        NSpatial: integer size of spatial dimension.
        NTemporal: integer size of temporal dimension.

    Returns:
        A list of length 4 with containing a square as possible hypercube.
    """

    restProc = numprocs;
    N = [0,0,0,0]
    for i in xrange(3):
        N[i] = NSpatial
        N[3] = NTemporal;
    while restProc >= 2:
        max_index = get_arg_max_index(N)
        N[max_index] /= 2
        restProc /= 2
        if (restProc < 2):
            break
    return N[::-1] # Reversing seems to be quicker.

def check_sub_dim_viability(subDims):
    """
    Checking if the sub dimensions are valid.

    Args:
        subDims: list of length 4 containing of the sub lattice dimensions.

    Raises:
        ValueError: exits if the length of the sub dimensions is not 4, or if 
            it is not containing only integers, or if the any sub dimension is 
            less than 2.
    """

    # Doing some basic error catching before the actual jobs starts
    if len(subDims) != 4 and sum([type(i) == int for i in subDims]) != 4:
        raise ValueError("%s is not a valid set of sub dimensions." % str(subDims))

    # Ensuring all dimensions are larger than 2
    for dim in subDims:
        if dim <= 2:
            raise ValueError("%d is not a valid dimension" % dim)

def set_field_configs(config, config_folder, config_start_number, base_folder=""):
    """
    Populates list of sorted field configs from folder config_folder into the 
    config dictionary.

    Args:
        config: dictionary containing the job setup.
        config_folder: string containing the inputFolder for where the 
            configurations will the loaded from.
        config_start_number: integer on which the the new configuration
            will be numbered from.

    Return:
        config: dictionary, now with the field config paths and start number.
    """

    config["load_field_configs"] = True
    config["inputFolder"] = os.path.normpath(config_folder)
    config["inputFolder"] = os.path.relpath(config["inputFolder"], base_folder)

    if os.path.isdir(config_folder):
        config["field_configs"] = natural_sort([fpath for fpath in os.listdir(config_folder) if (os.path.splitext(fpath)[-1] == ".bin")])
    else:
        raise OSError("Error: %s is not a directory." % (config_folder))

    corrected_configs = [] 
    for cfg in config["field_configs"]:
        cfg_number = [int(num) for num in re.split('(\d+)', cfg) if num.isdigit()][-1]
        if int(cfg_number) >= config_start_number:
            corrected_configs.append(cfg)
    config["field_configs"] = corrected_configs
    return config

def natural_sort(l):
    """
    Natural sorting function.

    Args:
        l: list of strings where each string contains a number, either on 
            format of 1,2,3... or 00001, 00002 ect.

    Returns:
        A sorted list.
    """

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)', key)]
    return sorted(l, key=alphanum_key)

class JobCreator:
    """Class for initializing jobs."""

    def __init__(self, dryrun, verbose):
        """
        Args:
            dryrun: bool, no folder or jobs will be initiated.
            verbose: a more verbose output.
        """

        self.dryrun = dryrun
        self.verbose = verbose
        self.CURRENT_PATH = os.getcwd()

        # Checking the .ids.json file for previous jobs and storing them in job-list.
        self.idFilesName = ".ids.json"
        if os.path.isfile(self.idFilesName):
            self.jobs = json.load(open(self.idFilesName, "r"))
        else:
            self.jobs = {}

    def _createDictionary(self, **kwargs):
        return_dict = {}
        ordered_dict_list = []
        for key, value in kwargs.items():
            return_dict[key] = value
        return return_dict

    def _create_folders(self):
        """Sets up the relevant folders for the job."""

        # Checking that we have an output folder.
        self._checkFolderPath(self.outputFolder)
        self._checkFolderPath(os.path.join(self.outputFolder, self.runName))
        if not self.uTest:

            if self.NFlows != 0:

                # Only in the case we are not creating any scalar fields
                if not self.create_fields_folders:
                    self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "flow_observables"))
                if self.create_fields_folders:
                    self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "scalar_fields"))

                for fobs in self.flow_observables:

                    # Only in the case we are not creating any scalar fields
                    if not self.create_fields_folders:
                        self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "flow_observables", fobs))

                    # The name of the observable for the scalar fields is different than the default names.
                    if self.create_fields_folders and fobs in AVAILABLE_SCALAR_FIELDS:
                        if "topc" in fobs.lower():
                            self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "scalar_fields", "topc"))
                        if "energy" in fobs.lower():
                            self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "scalar_fields", "energy"))

        if not self.load_field_configs and not self.create_fields_folders:
            self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "field_configurations"))
            if not self.uTest:
                self._checkFolderPath(os.path.join(self.outputFolder, self.runName, "observables"))

        self._checkFolderPath(os.path.join(self.inputFolder))
        self._checkFolderPath(os.path.join("input", self.runName))


    def _checkFolderPath(self, folder):
        """
        Function for checking if a folder exists, and if not creates on
        (unless we are doing a dryrun). Will be more verbose if prompted
        on initialization of class.

        Args:
            folder: folder path string.
        """

        if not os.path.isdir(os.path.join(self.base_folder, folder)):
            if not self.dryrun:
                os.mkdir(os.path.join(self.base_folder, folder))
            if self.dryrun or self.verbose:
                print "> mkdir %s" % os.path.join(self.base_folder, folder)

    def _check_lcfgr_path(self, config_dict):
        """Performs a check when loading and running from a configuration, 
        assuring the binary file we are loading exists at assumed location.
        """
        input_folder =  self._clean_file_path(config_dict["inputFolder"], quiet=True)
        if len(config_dict["load_config_and_run"]) != 0:
            cfg_path = os.path.join(os.path.normpath(self.base_folder) + 
                os.path.normpath(input_folder),
                config_dict["load_config_and_run"])
            if not os.path.isfile(cfg_path):
                exit("ERROR: binary file for 'load_config_and_run' not "
                    "found: {0:s}".format(cfg_path))


    def _clean_file_path(self, p, quiet=False):
        """
        Cleans file-path of base folder and of extraneous slashes.

        Args:
            p: file path string.
        Return:
            cleaned_p: cleaned file path string.
        """
        if p.startswith(os.path.abspath(self.base_folder) + os.sep):
            cleaned_p = os.path.normpath(os.sep + p + os.sep).replace("//", os.sep)
            cleaned_p = os.path.relpath(cleaned_p, self.base_folder)
        else:
            cleaned_p = p

        cleaned_p = (os.sep + os.path.normpath(cleaned_p) + os.sep).replace("//", os.sep)

        if self.verbose or self.dryrun and not quiet:
            print "Cleaned path from:\n    %s\nto\n    %s" % (p, cleaned_p)

        return cleaned_p

    def _create_json(self, config_dict):
        """
        Function that creates a json file for submitting to the c++ file.

        Args:
            config_dict: dictionary containing the configuration.
        """

        self.json_file_name = "config_%s.json" % self.runName
        json_dict = {}

        # Lattice related run variables
        json_dict["NSpatial"] = config_dict["N"]
        json_dict["NTemporal"] = config_dict["NT"]
        json_dict["subDims"] = config_dict["subDims"] # Will have to check for this being false
        json_dict["beta"] = config_dict["beta"]
        json_dict["NCf"] = config_dict["NCf"]
        json_dict["NCor"] = config_dict["NCor"]
        json_dict["NTherm"] = config_dict["NTherm"]
        json_dict["NFlows"] = config_dict["NFlows"]
        json_dict["NUpdates"] = config_dict["NUpdates"]

        # Data storage related variables
        json_dict["outputFolder"] = self._clean_file_path(config_dict["outputFolder"])
        json_dict["inputFolder"] = self._clean_file_path(config_dict["inputFolder"])
        json_dict["storeConfigurations"] = config_dict["storeCfgs"]
        json_dict["storeThermalizationObservables"] = config_dict["storeThermCfgs"]

        json_dict["load_config_and_run"] = config_dict["load_config_and_run"]
        json_dict["config_start_number"] = config_dict["config_start_number"]

        # Human readable output related variables
        json_dict["verbose"] = config_dict["verboseRun"]

        # Setup related variables
        json_dict["pwd"] = os.path.normpath(config_dict["base_folder"])
        json_dict["batchName"] = config_dict["runName"]
        json_dict["hotStart"] = config_dict["hotStart"]
        json_dict["RSTHotStart"] = config_dict["RSTHotStart"]
        json_dict["expFunc"] = config_dict["expFunc"]
        json_dict["action"] = config_dict["action"]
        json_dict["observables"] = config_dict["observables"]
        json_dict["flowObservables"] = config_dict["flowObservables"]
        json_dict["load_field_configs"] = config_dict["load_field_configs"]
        json_dict["chroma_config"] = config_dict["chroma_config"]
        json_dict["field_configs"] = config_dict["field_configs"]

        # Testing related variables
        json_dict["unitTesting"] = config_dict["uTest"]
        json_dict["unitTestingVerbose"] = config_dict["uTestVerbose"]
        json_dict["uTestFieldGaugeInvarince"] = config_dict["uTestFieldGaugeInvarince"]

        # Performance testing variables
        json_dict["performanceTesting"] = config_dict["performanceTesting"]
        json_dict["NExpTests"] = config_dict["NExpTests"]
        json_dict["NRandTests"] = config_dict["NRandTests"]
        json_dict["NDerivativeTests"] = config_dict["NDerivativeTests"]
        json_dict["TaylorPolDegree"] = config_dict["TaylorPolDegree"]

        # Data generation related variables
        json_dict["SU3Eps"] = config_dict["SU3Eps"]
        json_dict["flowEpsilon"] = config_dict["flowEpsilon"]
        json_dict["metropolisSeed"] = config_dict["metropolisSeed"]
        json_dict["randomMatrixSeed"] = config_dict["randomMatrixSeed"]
        json_dict["samplingFrequency"] = config_dict["samplingFrequency"]

        # Debugger
        json_dict["debug"] = config_dict["debug"]

        # Source
        json_fpath = os.path.join(self.base_folder, "input", self.json_file_name)

        # Dest
        json_cfgpath = os.path.join(self.base_folder, "input", self.runName, self.json_file_name)

        # Prints configuration file content if verbose or dryrun is true
        if self.dryrun or self.verbose:
            print "Writing json configuration file at location {0:<s}:\n".format(json_fpath)
            print json.dumps(json_dict, indent=4, separators=(", ", ": ")), "\n"

        # Creates configuration file
        if not self.dryrun:
            with file(json_fpath, "w+") as json_file:
                json.dump(json_dict, json_file, indent=4)
            shutil.copy(json_fpath, "%s.bak" % json_cfgpath)

    def submit_job(self, job_config, system, partition, excluded_nodes=False, ignore_tasks_per_node=False):
        if excluded_nodes:
            sbatch_exclusions = "#SBATCH --exclude=%s" % excluded_nodes
        else:
            sbatch_exclusions = ""

        # Checks if flow is sampling more observables than the regular config sampler, then sets it equal
        if (job_config["NFlows"] != 0):
            job_config["observables"] = job_config["flowObservables"]

        # Since topct gives us all of the observables basically for free,
        # we set the fobs to plaq, topc, energy, topct.
        if "topct" in job_config["flowObservables"]:
            job_config["flowObservables"] += ["plaq", "topc", "energy", "topct"]
            job_config["flowObservables"] = list(set(job_config["flowObservables"]))

        if "topct" in job_config["observables"]:
            job_config["observables"] += ["plaq", "topc", "energy", "topct"]
            job_config["observables"] = list(set(job_config["observables"]))

        # Retrieving config contents
        self.base_folder        = job_config["base_folder"]
        binary_filename         = job_config["bin_fn"]
        self.runName            = job_config["runName"]
        self.load_field_configs = job_config["load_field_configs"]
        threads                 = job_config["threads"]
        beta                    = job_config["beta"]
        NSpatial                = job_config["N"]
        NTemporal               = job_config["NT"]
        NTherm                  = job_config["NTherm"]
        NCor                    = job_config["NCor"] 
        NCf                     = job_config["NCf"]
        self.NFlows             = job_config["NFlows"]
        NUpdates                = job_config["NUpdates"]
        SU3Eps                  = job_config["SU3Eps"]
        self.flow_observables   = job_config["flowObservables"]
        self.inputFolder        = job_config["inputFolder"]
        self.outputFolder       = job_config["outputFolder"]
        self.observables        = job_config["observables"]
        flowEpsilon             = job_config["flowEpsilon"]
        storeCfgs               = job_config["storeCfgs"]
        storeThermCfgs          = job_config["storeThermCfgs"]
        hotStart                = job_config["hotStart"]
        RSTHotStart             = job_config["RSTHotStart"]
        subDims                 = job_config["subDims"]
        verboseRun              = job_config["verboseRun"]
        self.uTest              = job_config["uTest"]
        uTestVerbose            = job_config["uTestVerbose"]
        cpu_approx_runtime_hr   = job_config["cpu_approx_runtime_hr"]
        cpu_approx_runtime_min  = job_config["cpu_approx_runtime_min"]
        self.create_fields_folders = job_config["scalar_fields_folders"]
        self.user_mail          = job_config["user_mail"]

        self._check_lcfgr_path(job_config)

        # Checks that binary file exists in expected location
        if not os.path.isfile(os.path.join(self.CURRENT_PATH, binary_filename)):
            exit("ERROR: binary file path not in expected location %s/%s" % (
                self.CURRENT_PATH, binary_filename))

        # Ensures that we have a viable number of sub dimensions
        if len(subDims) != 0:
            check_sub_dim_viability(subDims)

        # Creates relevant folders
        self._create_folders()

        # Creates json config file
        self._create_json(job_config)

        # If we are on local computer(e.g. laptop), will create configuration file and quit
        if system == "local":
            sys.exit("Configuration file %s for local production created." 
                % os.path.join(self.base_folder, "input", self.json_file_name))

        # Setting job name before creating content file.
        job_name = "b{0:<3.2f}_cfg{1:<d}_nf{2:<d}_{3:<d}cube{4:<d}_{5:<d}threads".format(
            beta, NCf, self.NFlows, NSpatial, NTemporal, threads)

        # Setting approximated run time
        estimated_time = "{0:0>2d}:{1:0>2d}:00".format(
            cpu_approx_runtime_hr, cpu_approx_runtime_min)

        # # Choosing system
        # if system == "smaug":
        #     # Smaug batch file.
        #     # #SBATCH --exclude=smaug-b2 excludes an unstable node
        #     content = "#!/bin/bash"
        #     content += "\n#SBATCH --job-name={0:<s}".format(job_name)
        #     content += "\n#SBATCH --partition={0:<s}".format(partition)
        #     content += "\n#SBATCH --ntasks={0:<d}".format(threads)
        #     content += "\n#SBATCH --time={0:<s}".format(estimated_time)
        #     content += "\n" + sbatch_exclusions

        if system == "slurm":
            # Abel specific commands
            cpu_memory = job_config["cpu_memory"]
            account_name = job_config["account_name"]
            tasks_per_node = 16 # Maximum number of threads per node
            
            nodes = 1
            if threads > tasks_per_node:
                nodes = threads/tasks_per_node
            if threads % tasks_per_node != 0:
                if ignore_tasks_per_node:
                    print "Warning: Tasks(number of threads) are not divisible by 16: {0:d} % {1:d} = {2:d}".format(
                        threads, tasks_per_node, threads % tasks_per_node)
                else:
                    raise ValueError("Tasks(number of threads) have to be divisible by 16.")

            content = "#!/bin/bash"
            content += "\n#SBATCH --job-name={0:<s}".format(job_name)
            content += "\n#SBATCH --account={0:<s}".format(account_name)
            content += "\n#SBATCH --time={0:<s}".format(estimated_time)
            content += "\n#SBATCH --mem-per-cpu={0:<4d}".format(cpu_memory)
            content += "\n#SBATCH --partition={0:<s}".format(partition)
            content += "\n#SBATCH --ntasks={0:<d}".format(threads)
            content += "\n#SBATCH --nodes={0:<1d}".format(nodes)
            content += "\n#SBATCH --mail-type=BEGIN,TIME_LIMIT_10,END"
            content += "\n" + sbatch_exclusions + "\n"
            content += "\nsource /cluster/bin/jobsetup\n"
            content += "\nmodule purge                # clear any inherited modules"
            content += "\nmodule load openmpi.gnu     # loads mpi\n"
            content += "\nset -o errexit              # exit on errors\n"

        elif system == "torque":
            account_name = "ptg"
            tasks_per_node = 28
            
            nodes = 1
            if threads > tasks_per_node:
                nodes = threads/tasks_per_node

            if threads % (tasks_per_node * nodes) != 0:
                nodes += 1
            
            threads = tasks_per_node * nodes

            cpu_memory = job_config["cpu_memory"]
            # sys.exit("Error: not implemented setup for system: %s" % system) # Also, has to load module load Qt/5.6.2

            content = "#! /bin/bash -login"
            content += "\n#PBS -A {0:<s}".format(account_name)
            content += "\n#PBS -l walltime={0:<s},nodes={1:<1d}:ppn={2:<d},mem={3:<4d}MB".format(estimated_time, nodes, tasks_per_node, cpu_memory*threads)
            content += "\n#PBS -N {0:<s}".format(job_name)
            content += "\n#PBS -M {0:<s}".format(self.user_mail)
            content += "\n#PBS -m bea"
            # content += "\nmodule load GNU/6.2"
            # content += "\nmodule load OpenMPI/2.0.2"
            content += "\nmodule load GNU/4.9"
            content += "\nmodule load OpenMPI/1.10.0"
            #content += "\nmodule load Qt/5.6.2"

        elif system == "local":
            sys.exit("ERROR: this is a local production run. Should never see this error message.")
        else:
            sys.exit("ERROR: system %s not recognized." % system)

        # Setting run-command. Not needed anymore
        # if not system == "laconia":
        #     run_command = "mpirun -n {0:<d}".format(threads)
        # else :
        #     run_command = "mpirun -np {0:<d}".format(threads)

        run_command = "mpirun -np {0:<d}".format(threads)
        run_command += " "
        run_command += os.path.join(self.CURRENT_PATH, binary_filename)
        run_command += " "
        run_command += os.path.join(self.base_folder, "input", self.json_file_name)

        content += "\n" + run_command

        if system == "slurm":
            job = 'jobfile.slurm'
        else:
            job = 'jobfile.qsub'

        # Writes slurm file to be submitted
        if not self.dryrun:
            outfile = open(job, 'w')
            outfile.write(content)
            outfile.close()
        if self.dryrun or self.verbose:
            print "Writing %s to slurm batch file: \n\n" % job, content, "\n"

        # Sets up command based on system we have.
        if system == "slurm":
            cmd = ['sbatch', os.path.join(self.CURRENT_PATH,job)]
        else:
            cmd = ['qsub', os.path.join(self.CURRENT_PATH,job)]

        # Submits job
        if self.dryrun or self.verbose:
            print "> %s %s" % tuple(cmd)
            ID = 123456
        if not self.dryrun:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
            # tmp = "123 456" # Used for bug-hunting
            try:
                if system == "slurm":
                    ID = int(tmp.split()[-1]) # ID of job
                else:
                    tmp2 = tmp.split()[-1]
                    tmp2 = tmp2.split(".")[0]
                    ID = int(tmp2) # ID of laconia job
            except IndexError:
                print "ERROR: IndexError for line: \n", tmp, "--> exiting", exit(0)

        # Stores job in job dictionary
        job_dict =  self._createDictionary( 
            Partition = [0, 9, partition],
            RunName = [1, 35, self.runName],
            Beta = [2, 5, beta],
            N = [3, 4, NSpatial],
            NT = [4, 4, NTemporal],
            NCf = [5, 4, NCf],
            NTherm = [6, 6, NTherm],
            NCor = [7, 4, NCor],
            NUpdates = [8, 9, NUpdates],
            NFlows = [9, 7, self.NFlows],
            SU3Eps = [10, 6, SU3Eps],
            Threads = [11, 8, threads],
            StoreCfgs = [12, 10, bool(storeCfgs)],
            StoreThermCfgs = [13, len("StoreThermCfgs") + 1, bool(storeThermCfgs)],
            HotStart = [14, len("HotStart") + 1, bool(hotStart)],
            SubDims = [15, 15, ' '.join(map(str, subDims))],
            CPU_hr = [16, 7, cpu_approx_runtime_hr],
            CPU_min = [17, 7, cpu_approx_runtime_min])

        if not self.dryrun:
            self.jobs[ID] = job_dict

        # Changes name of job script
        if not self.dryrun:
            os.system('mv %s job_%d.sh' % (job, ID))
        if self.dryrun or self.verbose:
            print '> mv %s job_%d.sh' % (job, ID)

        # Updates ID file only if it is not a dryrun
        if not self.dryrun:
            self.update_id_file()

    def update_id_file(self):
        """Rewrites the .ids.txt file with additional job ID."""

        if self.dryrun:
            print "Updating %s file." % self.idFilesName
        else:
            if self.verbose:
                print "Updating %s file." % self.idFilesName
            with open(self.idFilesName, "w+") as f:
                json.dump(self.jobs, f, indent=4)

    def __select_system_cancel(self, system):
        """Selects system method for canceling HPC job."""
        if system == "torque":
            return "qdel"
        elif system == "slurm":
            return "scancel"
        else:
            raise KeyError("{} not a recognized system".format(system))


    def cancel_job(self, jobID, system):
        """Cancels jobID."""
        system_cmd = self.__select_system_cancel(system)
        if not self.dryrun:
            os.system("{0:s} {1:d}".format(system_cmd, jobID))
        if self.dryrun or self.verbose:
            print "> {0:s} {1:d}".format(system_cmd, jobID)

    def cancel_all_jobs(self):
        """Cancels all jobs."""
        system_cmd = self.__select_system_cancel(system)
        for i in self.jobs:
            if not self.dryrun:
                os.system("{0:s} {1:d}".format(system_cmd, int(i)))
            if self.dryrun or self.verbose:
                print "> {0:s} {1:d}".format(system_cmd, int(i))

    def list_jobs(self):
        """Shows jobs with ID numbers."""

        if len(self.jobs) == 0:
            print "No jobs running"
        else:
            # Takes the jobs out of their dictionary, and zips values and keys together for creating a header
            # Sorts jobs by their ID. Newest are printed at the bottom
            sorted_job_keys = sorted(self.jobs, key=lambda i: int(i))

            # Will order list based on the last header element, e.g. the most updated one
            last_job = self.jobs[sorted_job_keys[-1]]
            last_job_sorted = sorted(zip(last_job.keys(), last_job.values()), key=lambda i: i[-1][0])

            # Creates the header
            print "{0:<{w}}".format("ID", w=10),

            # Sorts based on the number give nto the header element
            for elem in last_job_sorted:
                # Retrieves the name of the header column
                name = elem[0]

                # Retrieves the width of the header column
                width = elem[1][1]
                
                # Prints the column name
                print "{0:<{w}}".format(name, w=width),

            # Prints the jobs
            for jobID in sorted_job_keys:
                # Prints job ID
                print "\n{0:<{w}}".format(jobID, w=10),

                # Takes the jobs out of their dictionary, and zips values and keys together for printing their values
                for i,item in enumerate(sorted(zip(self.jobs[jobID].keys(), self.jobs[jobID].values()), key=lambda i: i[-1][0])):
                    # Width based on the last job committed
                    width = last_job_sorted[i][1][1]

                    # Prints the item for the job id
                    print "{0:<{w}}".format(item[1][-1], w=width),

    def print_job_id_info(self, jobID):
        """
        Prints info for given job.

        Args:
            jobID: integer of the job id for a given system.
        """
        jobID = str(jobID)
        if len(self.jobs) == 0:
            print "No jobs running"
        else:
            # Takes the jobs out of their dictionary, and zips values and keys together for creating a header
            sorted_jobs = sorted(zip(self.jobs.values()[0].keys(), self.jobs.values()[0].values()), key=lambda i: i[-1][0])
            print "{0:<{w}}".format("ID", w=10),
            for i in sorted_jobs:
                print "{0:<{w}}".format(i[0], w=i[-1][1]),

            # Prints a single job
            print "\n{0:<{w}}".format(jobID, w=10),

            # Takes the jobs out of their dictionary, and zips values and keys together for printing their values
            for item in sorted(zip(self.jobs[jobID].keys(), self.jobs[jobID].values()), key=lambda i: i[-1][0]):
                print "{0:<{w}}".format(item[-1][-1], w=item[-1][1]),

    def clear_id_file(self):
        """Clears the ID file."""

        self.jobs = {}
        if not self.dryrun:
            self.update_id_file()
        if self.dryrun or self.verbose:
            print "Clearing ID file."

def main(args):
    # Default configuration file.
    if not os.path.isdir(os.path.join(os.getcwd(), "build")):
        raise EnvironmentError("Build folder is not present at location %s." % os.path.join(os.getcwd(), "build"))

    # Default config
    config_default = {
        "bin_fn"                    : "build-release/GLAC",
        "runName"                   : "defaultRun",
        "N"                         : 8, # Small lattice as default
        "NT"                        : 16,
        "subDims"                   : [],
        "beta"                      : 6.0,
        "NCf"                       : 100,
        "NCor"                      : 20,
        "NTherm"                    : 200,
        "NFlows"                    : 0,
        "NUpdates"                  : 10,
        "storeCfgs"                 : True,
        "storeThermCfgs"            : False,
        "verboseRun"                : False,
        "hotStart"                  : False,
        "RSTHotStart"               : False,
        "expFunc"                   : "morningstar",
        "action"                    : "wilson", # options: wilson, luscher - both first order
        "observables"               : ["plaq"],
        "flowObservables"           : ["plaq","topc","energy"],
        "load_field_configs"        : False,
        "load_config_and_run"       : "",
        "config_start_number"       : 0,
        "chroma_config"             : False,
        "base_folder"               : os.getcwd(),
        "inputFolder"               : "input",
        "outputFolder"              : "output",
        "field_configs"             : [],
        "uTest"                     : False,
        "uTestVerbose"              : False,
        "uTestFieldGaugeInvarince"  : "",
        "performanceTesting"        : False,
        "NExpTests"                 : int(1e6),
        "NRandTests"                : int(1e6),
        "NDerivativeTests"          : int(1e2),
        "TaylorPolDegree"           : 8,
        "SU3Eps"                    : 0.24,
        "flowEpsilon"               : 0.01,
        "metropolisSeed"            : 0,
        "randomMatrixSeed"          : 0,
        "threads"                   : 64,
        "scalar_fields_folders"     : False,
        "samplingFrequency"         : 25,
        "debug"                     : False,
        "cpu_approx_runtime_hr"     : 2, # In order to catch if config we are loading contains cpu approx time
        "cpu_approx_runtime_min"    : 0,
        "cpu_memory"                : 3800,
        "account_name"              : "nn2977k",
        "user_mail"                 : "h.m.m.vege@fys.uio.no",
    }

    ######## Initiating command line parser ########
    description_string = '''
    Program for starting large parallel Lattice Quantum Chromo Dynamics jobs.
    '''
    parser = argparse.ArgumentParser(prog='GLAC job creator', description=description_string)
    
    ######## Prints program version if prompted ########
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.2')
    parser.add_argument('--dryrun', default=False, action='store_true', help='Dryrun to not perform any critical actions.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='A more verbose output when generating.')

    ######## Sets up subparsers ########
    subparser = parser.add_subparsers(dest='subparser')

    ######## Job control ########
    sbatch_parser = subparser.add_parser('sbatch', help='Views, stops, clears and list jobs.')
    sbatch_group = sbatch_parser.add_mutually_exclusive_group(required=True)
    sbatch_group.add_argument('--scancel',                      default=False,      type=int, help='Cancel a job of given ID.')
    sbatch_group.add_argument('--scancel_all',                  default=False,      action='store_true', help='Cancel all jobs')
    sbatch_group.add_argument('-ls', '--list_jobs',             default=False,      action='store_true', help='List all jobs currently running.')
    sbatch_group.add_argument('--clear_id_file',                default=False,      action='store_true', help='Clears the job ID file.')
    sbatch_group.add_argument('-id', '--list_job_id',           default=False,      type=int, help='Shows details about job with given ID.')

    ######## Manual job setup ########
    job_parser = subparser.add_parser('setup', help='Sets up the job.')
    job_parser.add_argument('system',                           default=False,                                      type=str, choices=AVAILABLE_HPC_SYSTEMS, help='Specify system we are running on.')
    job_parser.add_argument('threads',                          default=False,                                      type=int, help='Number of threads to run on')
    job_parser.add_argument('-p', '--partition',                default="normal",                                   type=str, help='Specify partition to run program on.')
    job_parser.add_argument('-rn', '--run_name',                default=config_default["runName"],                  type=str, help='Specify the run name')

    # Lattice related run variables
    job_parser.add_argument('-N', '--NSpatial',                 default=config_default["N"],                        type=int, help='spatial lattice dimension')
    job_parser.add_argument('-NT', '--NTemporal',               default=config_default["NT"],                       type=int, help='temporal lattice dimension')
    job_parser.add_argument('-sd', '--subDims',                 default=False,                                      type=int, nargs=4, help='List of sub lattice dimension sizes, length 4')
    job_parser.add_argument('-b', '--beta',                     default=config_default["beta"],                     type=float, help='beta value')
    job_parser.add_argument('-NCfgs', '-Ncfg', '-NCf', '--NConfigs', default=config_default["NCf"],                      type=int, help='number of configurations to generate')
    job_parser.add_argument('-NCor', '-NCorr', '--NCor',        default=config_default["NCor"],                     type=int, help='number of correlation updates to perform')
    job_parser.add_argument('-NTh', '--NTherm',                 default=config_default["NTherm"],                   type=int, help='number of thermalization steps')
    job_parser.add_argument('-NFlows', '--NFlows',              default=config_default["NFlows"],                   type=int, help='number of flows to perform per configuration')
    job_parser.add_argument('-NUp', '--NUpdates',               default=config_default["NUpdates"],                 type=int, help='number of updates per link')

    # Data storage related variables
    job_parser.add_argument('-sc', '--storeCfgs',               default=config_default["storeCfgs"],                type=int, choices=[0,1], help='Specifying if we are to store configurations')
    job_parser.add_argument('-st', '--storeThermCfgs',          default=config_default["storeThermCfgs"],           type=int, choices=[0,1], help='Specifies if we are to store the thermalization plaquettes')
    job_parser.add_argument('-bf', '--base_folder',             default=config_default["base_folder"],              type=str, help='Sets the base folder. Default is os.path.getcwd().') # Human readable output related variables
    job_parser.add_argument('-vr', '--verboseRun',              default=config_default["verboseRun"],               action='store_true', help='Verbose run of GLAC. By default, it is off.')

    # Setup related variables
    job_parser.add_argument('-hs', '--hotStart',                default=config_default["hotStart"],                 type=int, choices=[0,1], help='Hot start or cold start')
    job_parser.add_argument('-rsths', '--RSTHotStart',          default=config_default["RSTHotStart"],              type=int, choices=[0,1], help='RST hot start is closer to unity')
    job_parser.add_argument('-expf', '--expFunc',               default=config_default["expFunc"],                  type=str, choices=AVAILABLE_EXP_FUNCS, help='Sets the exponentiation function to be used in flow. Default is method by Morningstar.')
    job_parser.add_argument('-obs', '--observables',            default=config_default["observables"],              type=str, choices=AVAILABLE_OBSERVABLES, nargs='+', help='Observables to sample for in flow.')
    job_parser.add_argument('-fobs', '--flowObservables',       default=config_default["flowObservables"],          type=str, choices=AVAILABLE_OBSERVABLES, nargs='+', help='Observables to sample for in flow.')
    job_parser.add_argument('-act', '--action',                 default=config_default["action"],                   type=str, choices=AVAILABLE_ACTIONS, nargs=None, help='Action to use.')

    # Data generation related variables
    job_parser.add_argument('-SU3Eps', '--SU3Epsilon',          default=config_default["SU3Eps"],                   type=float, help='SU3 epsilon random increment value.')
    job_parser.add_argument('-fEps', '--flowEpsilon',           default=config_default["flowEpsilon"],              type=float, help='Flow epsilon derivative small change value.')
    job_parser.add_argument('-mSeed', '--metropolisSeed',       default=config_default["metropolisSeed"],           type=float, help='Seed for the Metropolis algorithm.')
    job_parser.add_argument('-rSeed', '--randomMatrixSeed',     default=config_default["randomMatrixSeed"],         type=float, help='Seed for the random matrix generation.')

    # Other usefull parsing options
    job_parser.add_argument('-sq', '--square',                  default=False,                                      action='store_true', help='Enforce square sub lattices(or as close as possible).')
    job_parser.add_argument('-chr', '--cpu_approx_runtime_hr',  default=2,                                          type=int, help='Approximate cpu time in hours that will be used')
    job_parser.add_argument('-cmin', '--cpu_approx_runtime_min',default=0,                                          type=int, help='Approximate cpu time in minutes that will be used')
    job_parser.add_argument('-ex', '--exclude',                 default=False,                                      type=str, nargs='+', help='Nodes to exclude.')
    job_parser.add_argument('-lcfg', '--load_configurations',   default=config_default["load_field_configs"],       type=str, help='Loads configurations from a folder by scanning and for files with .bin extensions.')
    job_parser.add_argument('-chroma', '--chroma_config',       default=config_default["chroma_config"],            action='store_true', help='If flagged, loads the configuration as a chroma configuration.')
    job_parser.add_argument('-lcfgr', '--load_config_and_run',  default=False,                                      type=str, help='Loads a configuration that is already thermalized and continues generating N configurations based on required -NCfgs argument.')
    job_parser.add_argument('-cfgnum', '--config_start_number', default=config_default["config_start_number"],      type=int, help='Starts naming the configuration from this number.')
    job_parser.add_argument('-igntsk', '--ignore_tasks_per_node', default=False,                                    action='store_true', help='If enabled, will ignore requirement of having 16 tasks per node.')

    # Debug option
    job_parser.add_argument('--debug',                          default=config_default["debug"],                    action='store_true', help='Debug option. Will check lattices for corruption and zeros.')

    ######## Abel specific commands ########
    job_parser.add_argument('--cpu_memory',                     default=config_default["cpu_memory"],               type=int, help='CPU memory to be allocated to each core')
    job_parser.add_argument('--account_name',                   default=config_default["account_name"],             type=str, help='Account name associated to the abel cluster')

    ######## Job load parser ########
    load_parser = subparser.add_parser('load', help='Loads a configuration file into the program')
    load_parser.add_argument('file',                            default=False,                                      type=str, help='Loads config file')
    load_parser.add_argument('-s', '--system',                  default=False,                                      type=str, required=True, choices=AVAILABLE_HPC_SYSTEMS, help='Cluster name')
    load_parser.add_argument('-p', '--partition',               default="normal",                                   type=str, help='Partition to run on. Default is normal. If some nodes are down, manual input may be needed.')
    load_parser.add_argument('-lcfg', '--load_configurations',  default=config_default["load_field_configs"],       type=str, help='Loads configurations from a folder by scanning and for files with .bin extensions.')
    load_parser.add_argument('-lcfgr', '--load_config_and_run', default=False,                                      type=str, help='Loads a configuration that is already thermalized and continues generating N configurations based on required -NCfgs argument.')
    load_parser.add_argument('-NCfgs', '-NCfg', '-NCf', '--NConfigs', default=False,                                type=int, help='N configurations to generate based on loaded configuration.')
    load_parser.add_argument('-NFlows', '--NFlows',             default=False,                                      type=int, help='number of flows to perform per configuration')
    load_parser.add_argument('-chroma', '--chroma_config',      default=config_default["chroma_config"],            action='store_true', help='If flagged, loads the configuration as a chroma configuration.')
    load_parser.add_argument('-lhr', '--load_config_hr_time_estimate', default=None,                                type=int, help='Number of hours that we estimate we need to run the loaded configurations for.')
    load_parser.add_argument('-lmin', '--load_config_min_time_estimate', default=None,                              type=int, help='Approximate cpu time in minutes that will be used.')
    load_parser.add_argument('-bf', '--base_folder',            default=config_default["base_folder"],              type=str, help='Sets the base folder. Default is os.path.getcwd().')
    load_parser.add_argument('-nf', '--no_flow',                default=False,                                      action='store_true', help='If toggled, will not perform any flows.')
    load_parser.add_argument('-cfgnum', '--config_start_number',default=config_default["config_start_number"],      type=int, help='Starts naming the configuration from this number.')
    load_parser.add_argument('-rn', '--run_name',               default=False,                                      type=str, help='Specify the run name')
    load_parser.add_argument('-ex', '--exclude',                default=False,                                      type=str, nargs='+', help='Nodes to exclude.')
    load_parser.add_argument('-NUp', '--NUpdates',              default=False,                                      type=int, help='number of updates per link')
    load_parser.add_argument('-NCor', '-NCorr', '--NCor',       default=False,                                      type=int, help='number of correlation updates to perform')
    load_parser.add_argument('--debug',                         default=False,                                      action='store_true', help='Debug option. Will check lattices for corruption and zeros.')
    load_parser.add_argument('-vr', '--verboseRun',             default=config_default["verboseRun"],              action='store_true', help='Verbose run of GLAC. By default, it is off.')
    load_parser.add_argument('-igntsk', '--ignore_tasks_per_node', default=False,                                   action='store_true', help='If enabled, will ignore requirement of having 16 tasks per node.')
    load_parser.add_argument('--account_name',                   default=config_default["account_name"],             type=str, help='Account name associated to the abel cluster')

    ######## Unit test parser ########
    unit_test_parser = subparser.add_parser('utest', help='Runs unit tests embedded in the GLAC program. Will exit when complete.')
    unit_test_parser.add_argument('system',                     default=False,                                      type=str, choices=AVAILABLE_HPC_SYSTEMS, help='Specify system we are running on.')
    unit_test_parser.add_argument('threads',                    default=False,                                      type=int, help='Number of threads to run on')
    unit_test_parser.add_argument('-vr', '--verboseRun',        default=config_default["verboseRun"],               action='store_true', help='Prints more information during testing.')
    unit_test_parser.add_argument('-cgi', '--check_gauge_invariance', default=False,                                type=str, help='Loads and checks the gauge field invariance of a field.')
    unit_test_parser.add_argument('-sq', '--square',            default=False,                                      action='store_true', help='Enforce square sub lattices(or as close as possible).')
    unit_test_parser.add_argument('-N', '--NSpatial',           default=config_default["N"],                        type=int, help='spatial lattice dimension')
    unit_test_parser.add_argument('-NT', '--NTemporal',         default=config_default["NT"],                       type=int, help='temporal lattice dimension')
    unit_test_parser.add_argument('-ex', '--exclude',           default=False,                                      type=str, nargs='+', help='Nodes to exclude.')
    unit_test_parser.add_argument('-igntsk', '--ignore_tasks_per_node', default=False,                              action='store_true', help='If enabled, will ignore requirement of having 16 tasks per node.')

    ######## Performance test parser ########
    performance_test_parser = subparser.add_parser('perf_test', help='Runs performance tests on the certain components of the GLAC program. Will exit when complete.')
    performance_test_parser.add_argument('system',              default=False,                                      type=str, choices=AVAILABLE_HPC_SYSTEMS, help='Specify system we are running on.')
    performance_test_parser.add_argument('threads',             default=False,                                      type=int, help='Number of threads to run on')
    performance_test_parser.add_argument('-NExpTests',          default=config_default["NExpTests"],                type=int, help='Number of exponentiation tests we will run.')
    performance_test_parser.add_argument('-NRandTests',         default=config_default["NRandTests"],               type=int, help='Number of random tests we will run.')
    performance_test_parser.add_argument('-NDerivativeTests',   default=config_default["NDerivativeTests"],         type=int, help='Number of full lattice derivative tests we will run.')
    performance_test_parser.add_argument('-TaylorPolDegree',    default=config_default["TaylorPolDegree"],          type=int, help='Degree of the Taylor polynomial for exponentiation(default is 8).')
    performance_test_parser.add_argument('-ex', '--exclude',    default=False,                                      type=str, nargs='+', help='Nodes to exclude.')
    performance_test_parser.add_argument('-igntsk', '--ignore_tasks_per_node', default=False,                       action='store_true', help='If enabled, will ignore requirement of having 16 tasks per node.')

    ######## Field densities sampler ########
    field_density_parser = subparser.add_parser("field_density", help="Will run a single config through the energyTopcFieldDensity observable and produce observables of the entire field from them")
    
    # Setup related variable
    field_density_parser.add_argument('-s', '--system',        default=False,                                      type=str, required=True, choices=AVAILABLE_HPC_SYSTEMS, help='Cluster name')
    field_density_parser.add_argument('-nt', '--threads',      default=config_default["threads"],                  type=int, required=True, help='Number of threads to run on')
    field_density_parser.add_argument('-cfgf', '--config_file', default=False,                                     type=str, help='Loads config file')
    field_density_parser.add_argument('-p', '--partition',     default="normal",                                   type=str, help='Partition to run on. Default is normal. If some nodes are down, manual input may be needed.')
    field_density_parser.add_argument('-lcfg', '--load_configuration', default=False,                              type=str, required=True, help='Loads a single configuration with .bin extension.')
    field_density_parser.add_argument('-rn', '--run_name',     default=config_default["runName"],                  type=str, help='Specify the run name')
    field_density_parser.add_argument('-lhr', '--load_config_hr_time_estimate', default=None,                      type=int, help='Number of hours that we estimate we need to run the loaded configurations for.')
    field_density_parser.add_argument('-lmin', '--load_config_min_time_estimate', default=None,                    type=int, help='Approximate cpu time in minutes that will be used.')
    field_density_parser.add_argument('-ex', '--exclude',      default=False,                                      type=str, nargs='+', help='Nodes to exclude.')
    field_density_parser.add_argument('-igntsk', '--ignore_tasks_per_node', default=False,                         action='store_true', help='If enabled, will ignore requirement of having 16 tasks per node.')
    
    # Lattice related run variables
    field_density_parser.add_argument('-b', '--beta',          default=False,                                      type=float, help='beta value')
    field_density_parser.add_argument('-N', '--NSpatial',      default=False,                                      type=int, help='spatial lattice dimension')
    field_density_parser.add_argument('-NT', '--NTemporal',    default=False,                                      type=int, help='temporal lattice dimension')
    field_density_parser.add_argument('-sd', '--subDims',      default=config_default["subDims"],                  type=int, nargs=4, help='List of sub lattice dimension sizes, length 4')
    field_density_parser.add_argument('-NFlows', '--NFlows',   default=config_default["NFlows"],                   type=int, help='number of flows to perform per configuration')
    field_density_parser.add_argument('-fEps', '--flowEpsilon', default=config_default["flowEpsilon"],             type=float, help='Flow epsilon derivative small change value.')
    field_density_parser.add_argument('-sq', '--square',       default=False,                                      action='store_true', help='Enforce square sub lattices(or as close as possible).')
    
    # Data storage related variables
    field_density_parser.add_argument('-bf', '--base_folder',  default=config_default["base_folder"],              type=str, help='Sets the base folder. Default is os.path.getcwd().')
    field_density_parser.add_argument('-sf', '--samplingFrequency', default=config_default["samplingFrequency"],   type=int, help='Sets the sampling frequency of the flow. Lattice to be written to file every given number.')
    field_density_parser.add_argument('-vr', '--verboseRun',   default=config_default["verboseRun"],               action='store_true', help='Verbose run of GLAC. By default, it is off.')
    
    # Debug variables
    field_density_parser.add_argument('--debug',               default=config_default["debug"],                    action='store_true', help='Debug option. Will check lattices for corruption and zeros.')

    args = parser.parse_args()

    # args = parser.parse_args(['--dryrun', 'field_density', '-s', 'local', '-nt', '8', '-lcfg', 'output/testRunLocal/field_configurations/testRunLocal_beta6.000000_spatial8_temporal16_threads8_config00000.bin', '-rn', 'lattice_field_density_test'])

    dryrun = args.dryrun
    verbose = args.verbose

    # Initiates JobCreator class for running jobs ect
    s = JobCreator(dryrun, verbose)

    # Loads one configuration
    if args.subparser == 'load':
        """
        COMMAND FOR LOADING SCRIPT JOBS AND GENERATING CONFIGS OR FLOWING
        """
        configuration = ast.literal_eval(open(args.file, "r").read())

        # Sets the base folder
        configuration["base_folder"] = args.base_folder

        # If we are to load and flow configurations
        if args.load_configurations:
            if args.load_config_and_run:
                # Error catching in case user is using load_config_and_run together with load_configurations.
                sys.exit("ERROR: can not load and run configurations(-lcfgr) together with load configurations(-lcfg).")
            
            # Requiring flow to be specified if we are loading configurations to flow
            if configuration["NFlows"] == 0 or args.no_flow:
                sys.exit("ERROR: when loading configuration for to flow, need to specifiy number of flows.")

            # Sets the number of configurations to run for to be zero just in case
            configuration["NCf"] = 0

            # Requiring an new estimate of the run time if we are flowing
            configuration = set_field_configs(configuration, args.load_configurations, args.config_start_number, base_folder=configuration["base_folder"])
            if args.load_config_min_time_estimate == None or args.load_config_hr_time_estimate == None:
                if "cpu_approx_runtime_hr" in configuration and "cpu_approx_runtime_min" in configuration:
                    args.load_config_hr_time_estimate = configuration["cpu_approx_runtime_hr"]
                    args.load_config_min_time_estimate = configuration["cpu_approx_runtime_min"]
                else:
                    sys.exit("ERROR: Need an estimate of the runtime for the flowing of configurations.")

        # For loading and running configurations
        if args.load_config_and_run:
            if not args.NConfigs and configuration["NCf"]:
                # Error catching, as we require to know how many addition configurations we wish to create from loaded configuration.
                sys.exit("ERROR: we require to know how many addition configurations we wish to create from loaded configuration(specified by -lcfgr).")
            _config_path, _config_file = os.path.split(args.load_config_and_run)
            assert not os.path.isdir(args.load_config_and_run), "%s is a folder and not a binary(.bin) file." % args.load_config_and_run
            assert _config_file.split(".")[-1] == "bin", "%s is not a binary(.bin) file." % _config_file
            configuration["load_config_and_run"] = _config_file
            configuration["inputFolder"] = os.path.normpath(_config_path)

            configuration["NCf"] = args.NConfigs
            configuration["NTherm"] = 0

        # Excludes certain nodes if arguments have been provided
        if args.exclude:
            excluded_nodes = ','.join(args.exclude)
        else:
            excluded_nodes = ""

        # Populate configuration with default values if certain keys are not present
        if args.run_name != False:
            configuration["runName"] = args.run_name
        if args.NConfigs != False:
            configuration["NCf"] = args.NConfigs
        if args.NFlows != False:
            configuration["NFlows"] = args.NFlows
        if args.NCor != False:
            configuration["NCor"] = args.NCor
        if args.NUpdates != False:
            configuration["NUpdates"] = args.NUpdates
        if args.load_config_min_time_estimate != None:
            configuration["cpu_approx_runtime_min"] = args.load_config_min_time_estimate
        if args.load_config_hr_time_estimate != None:
            configuration["cpu_approx_runtime_hr"] = args.load_config_hr_time_estimate
        if args.no_flow:
            configuration["NFlows"] = 0
            configuration["flowObservables"] = []
        if args.debug != False:
            configuration["debug"] = args.debug
        for key in config_default.keys():
            if not key in configuration:
                configuration[key] = config_default[key]

        configuration["account_name"] = args.account_name


        configuration["verboseRun"] = args.verboseRun
        configuration["chroma_config"] = args.chroma_config
        configuration["config_start_number"] = args.config_start_number

        # Submitting job
        s.submit_job(configuration, args.system, args.partition, excluded_nodes, ignore_tasks_per_node=args.ignore_tasks_per_node)

    elif args.subparser == 'setup':
        """
        COMMAND FOR SETTING UP JOBS WITHOUT ANY PREDEFINED SCRIPTS
        """
        config_default["runName"]                   = args.run_name
        config_default["threads"]                   = args.threads
        config_default["N"]                         = args.NSpatial
        config_default["NT"]                        = args.NTemporal
        config_default["beta"]                      = args.beta
        if args.NConfigs == 0:
            config_default["NCf"]                   = 1
        else:
            config_default["NCf"]                   = args.NConfigs
        config_default["NCor"]                      = args.NCor
        config_default["NTherm"]                    = args.NTherm
        config_default["NFlows"]                    = args.NFlows
        config_default["NUpdates"]                  = args.NUpdates
        config_default["storeCfgs"]                 = bool(args.storeCfgs)
        config_default["storeThermCfgs"]            = bool(args.storeThermCfgs)
        config_default["verboseRun"]                = args.verboseRun
        config_default["hotStart"]                  = bool(args.hotStart)
        config_default["RSTHotStart"]               = bool(args.RSTHotStart)
        config_default["expFunc"]                   = args.expFunc
        config_default["action"]                    = args.action
        config_default["observables"]               = args.observables
        config_default["flowObservables"]           = args.flowObservables
        config_default["SU3Eps"]                    = args.SU3Epsilon
        config_default["flowEpsilon"]               = args.flowEpsilon
        config_default["metropolisSeed"]            = args.metropolisSeed
        config_default["randomMatrixSeed"]          = args.randomMatrixSeed
        config_default["cpu_approx_runtime_hr"]     = args.cpu_approx_runtime_hr
        config_default["cpu_approx_runtime_min"]    = args.cpu_approx_runtime_min
        config_default["account_name"]              = args.account_name
        config_default["cpu_memory"]                = args.cpu_memory
        config_default["base_folder"]               = args.base_folder
        config_default["debug"]                     = args.debug

        # Non-trivial default values
        if args.subDims:
            check_sub_dim_viability(args.subDims)
            config_default["subDims"] = args.subDims
        if args.square:
            config_default["subDims"] = create_square(config_default["threads"], config_default["N"], config_default["NT"])

        # Excludes certain nodes if arguments have been provided
        if args.exclude:
            excluded_nodes = ','.join(args.exclude)
        else:
            excluded_nodes = ""

        if args.load_config_and_run != False:
            _config_path, _config_file = os.path.split(args.load_config_and_run)
            assert not os.path.isdir(args.load_config_and_run), "%s is a folder and not a binary(.bin) file." % args.load_config_and_run
            assert _config_file.split(".")[-1] == "bin", "%s is not a binary(.bin) file." % _config_file
            config_default["load_config_and_run"] = _config_file
            config_default["inputFolder"] = os.path.normpath(_config_path)

        config_default["config_start_number"] = args.config_start_number

        if args.load_configurations:
            config_default = set_field_configs(config_default, args.load_configurations, args.config_start_number,  base_folder=config_default["base_folder"])
            config_default["chroma_config"] = args.chroma_config

            # Requiring flow to be specified if we are loading configurations to flow
            if config_default["NFlows"] == 0:
                sys.exit("ERROR: when loading configuration for to flow, need to specifiy number of flows.")

        # Submitting job
        s.submit_job(config_default, args.system, args.partition, excluded_nodes, ignore_tasks_per_node=args.ignore_tasks_per_node)

    elif args.subparser == 'sbatch':
        """
        COMMAND FOR VIEWING JOBS
        """
        if args.scancel:
            s.cancel_job(args.scancel, args.system)
            return
        if args.scancel_all:
            s.cancel_all_jobs(args.system)
            return
        if args.list_jobs:
            s.list_jobs()
            return
        if args.clear_id_file:
            s.clear_id_file()
            return
        if args.list_job_id != False:
            s.print_job_id_info(args.list_job_id)
            return

    elif args.subparser == 'utest':
        """
        COMMAND FOR UNIT TESTING
        """
        config_default["runName"] = "defaultTestRun"
        config_default["threads"] = args.threads
        config_default["uTest"] = True
        config_default["cpu_approx_runtime_hr"] = 0
        config_default["cpu_approx_runtime_min"] = 20

        # In order to create output field cfg folders
        config_default["load_field_configs"] = False
        config_default["scalar_fields_folders"] = False

        partition = "normal"
        system = args.system
        config_default["uTestVerbose"] = args.verboseRun
        if args.check_gauge_invariance:
            config_default["uTestFieldGaugeInvarince"] = args.check_gauge_invariance
            if not args.NSpatial or not args.NTemporal:
                sys.exit("ERROR: need to specifiy dimensions of loaded lattice.")
        
        config_default["N"] = args.NSpatial
        config_default["NT"] = args.NTemporal
        if args.square:
            config_default["subDims"] = create_square(config_default["threads"], config_default["N"], config_default["NT"])

        # Checks if we are to exclude any of the nodes
        if args.exclude:
            excluded_nodes = ','.join(args.exclude)
        else:
            excluded_nodes = ""

        # Submitting job
        s.submit_job(config_default, system, partition, excluded_nodes, ignore_tasks_per_node=args.ignore_tasks_per_node)

    elif args.subparser == 'perf_test':
        """
        COMMAND FOR RUNNING PERFORMANCE TESTS.
        """
        config_default["runName"] = "defaultPerformanceRun"
        config_default["performanceTesting"] = True
        config_default["cpu_approx_runtime_hr"] = 0
        config_default["cpu_approx_runtime_min"] = 20
        config_default["N"] = 8
        config_default["NT"] = 16
        config_default["threads"] = args.threads
        config_default["subDims"] = create_square(config_default["threads"], config_default["N"], config_default["NT"])
        partition = "normal"
        config_default["NExpTests"] = args.NExpTests
        config_default["NRandTests"] = args.NRandTests
        config_default["NDerivativeTests"] = args.NDerivativeTests
        config_default["TaylorPolDegree"] = args.TaylorPolDegree

        # Checks if we are to exclude any of the nodes
        if args.exclude:
            excluded_nodes = ','.join(args.exclude)
        else:
            excluded_nodes = ""

        # Submitting job
        s.submit_job(config_default, args.system, partition, excluded_nodes, ignore_tasks_per_node=args.ignore_tasks_per_node)

    elif args.subparser == "field_density":
        field_density_parser = subparser.add_parser("field_density", help=("Will run a single config through the "
            "energyTopcFieldDensity observable and produce observables of the entire field from them"))

        if args.config_file != False:
            # Loads the config file
            configuration = ast.literal_eval(open(args.config_file, "r").read())

            # Ensures we do not have a conflicting job setup
            if args.NSpatial != configuration["N"] and args.NSpatial:
                raise KeyError("Spatial dimension already provided in configuration file %s." % args.config_file)
            if args.NTemporal != configuration["NT"] and args.NTemporal:
                raise KeyError("Temporal dimension already provided in configuration file %s." % args.config_file)
            if args.beta != configuration["beta"] and args.beta:
                raise KeyError("Beta value already provided in configuration file %s." % args.config_file)
            if args.subDims != configuration["subDims"] and args.subDims != config_default["subDims"]:
                raise KeyError("Sub-dimensions is already provided in configuration file %s." % args.config_file)
        else:
            # Sets default configuration as baseline
            configuration = config_default

            # If no config file is provided, we require N, NT and beta to be provided
            if not args.NSpatial or not args.NTemporal or not args.beta:
                raise KeyError("please provide dimensions N, NT and beta.")
            configuration["N"] = args.NSpatial
            configuration["NT"] = args.NTemporal
            configuration["beta"] = args.beta

            # Non-trivial default values
            if len(args.subDims):
                check_sub_dim_viability(args.subDims)
                config_default["subDims"] = args.subDims
            if args.square:
                config_default["subDims"] = create_square(config_default["threads"], config_default["N"], config_default["NT"])

        # Creates a folder to place scalar field values in.
        configuration["scalar_fields_folders"] = True

        # Sets the number of configurations to run for to be zero just in case
        configuration["NCf"] = 0

        configuration["debug"] = args.debug

        # Requiring an new estimate of the run time if we are flowing
        _config_path, _config_file = os.path.split(args.load_configuration)
        assert not os.path.isdir(args.load_configuration), "%s is a folder and not a binary(.bin) file." % args.load_configuration
        assert _config_file.split(".")[-1] == "bin", "%s is not a binary(.bin) file." % _config_file
        config_default["load_config_and_run"] = _config_file
        config_default["inputFolder"] = os.path.normpath(_config_path)
        config_default["config_start_number"] = 0

        # Sets some regular variables
        configuration["threads"] = args.threads
        configuration["runName"] = args.run_name
        configuration["NFlows"] = args.NFlows
        configuration["base_folder"] = args.base_folder
        configuration["flowEpsilon"] = args.flowEpsilon
        configuration["verboseRun"] = args.verboseRun
        configuration["observables"] = ["energyTopcFieldDensity"]
        configuration["flowObservables"] = ["energyTopcFieldDensity"]
        configuration["samplingFrequency"] = args.samplingFrequency

        # Sets nodes to skip, usefull when clusters contain bad nodes
        if args.exclude:
            excluded_nodes = ','.join(args.exclude)
        else:
            excluded_nodes = ""

        # Requires an estimation of the total run time when not running on local computer
        if args.system != "local":
            if args.load_config_min_time_estimate == None or args.load_config_hr_time_estimate == None:
                sys.exit("ERROR: Need an estimate of the runtime for the flowing of configurations.")
            configuration["cpu_approx_runtime_min"] = args.load_config_min_time_estimate
            configuration["cpu_approx_runtime_hr"] = args.load_config_hr_time_estimate

        for key in config_default.keys():
            if not key in configuration:
                configuration[key] = config_default[key]

        s.submit_job(configuration, args.system, args.partition, excluded_nodes, ignore_tasks_per_node=args.ignore_tasks_per_node)

    else:
        """
        ERROR CATCHING
        """
        print 'Parse error: %s \n--> exiting' % args
        exit(0)

if __name__ == '__main__':
    main(sys.argv[1:])