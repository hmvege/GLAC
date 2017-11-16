import os, subprocess, time, sys, argparse, json

def getArgMaxIndex(N):
    # For getting the maximum index of an list.
    val = N[0]
    index = 0
    for i in xrange(4):
        if N[i] > val:
            val = N[i]
            index = i
    return index

def createSquare(numprocs,NSpatial,NTemporal):
    # Create a square sub lattice
    restProc = numprocs;
    N = [0,0,0,0]
    for i in xrange(3):
        N[i] = NSpatial
        N[3] = NTemporal;
    while restProc >= 2:
        max_index = getArgMaxIndex(N)
        N[max_index] /= 2
        restProc /= 2
        if (restProc < 2):
            break
    return N[::-1] # Reversing seems to be quicker.

def checkSubDimViability(subDims):
    # Doing some basic error catching before the actual jobs starts
    if len(subDims) != 4 and sum([type(i) == int for i in subDims]) != 4:
        raise ValueError("%g is not a valid set of sub dimensions." % subDims)
    # Ensuring all dimensions are larger than 2
    for dim in subDims:
        if dim <= 2: exit("Error: %d is not a valid dimension" % dim)


class Slurm:
    def __init__(self, dryrun):
        self.dryrun = dryrun
        self.CURRENT_PATH = os.getcwd()
        # Checking that we have an output folder.
        self._checkFolderPath('output')
        self._checkFolderPath('output/flow_observables')
        self._checkFolderPath('output/field_configurations')
        self._checkFolderPath('output/observables')
        self._checkFolderPath('input')
        # Checking the .ids.json file for previous jobs and storing them in job-list.
        self.idFilesName = '.ids.json'
        if os.path.isfile(self.idFilesName) and not self.dryrun:
            self.jobs = json.load(open(self.idFilesName,"r"))
        else:
            self.jobs = {}

    def _checkFolderPath(self, folder):
        # Function for checking if a folder exists, and if not creates on(unless we are doing a dryrun)
        if not os.path.isdir('%s/%s' % (self.CURRENT_PATH,folder)):
            if self.dryrun:
                print '> mkdir %s/%s' % (self.CURRENT_PATH,folder)
            else:
                os.mkdir('%s/%s' % (self.CURRENT_PATH,folder))

    def _create_json(self,config_dict):
        # Function that creates a json file for submitting to the c++ file.
        with file("%s/input/config.json" % self.CURRENT_PATH,"w") as json_file:
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
            # Data storage related variables
            json_dict["storeConfigurations"] = config_dict["storeCfgs"]
            json_dict["storeThermalizationObservables"] = config_dict["storeThermCfgs"]
            # Human readable output related variables
            json_dict["verbose"] = config_dict["verboseRun"]
            # Setup related variables
            json_dict["pwd"] = self.CURRENT_PATH
            json_dict["batchName"] = config_dict["runName"]
            json_dict["hotStart"] = config_dict["hotStart"]
            json_dict["exponentiationFunction"] = config_dict["exponentiationFunction"]
            json_dict["observables"] = config_dict["observables"]
            json_dict["flowObservables"] = config_dict["flowObservables"]
            # Testing related variables
            json_dict["unitTesting"] = config_dict["uTest"]
            json_dict["unitTestingVerbose"] = config_dict["uTestVerbose"]
            # Data generation related variables
            json_dict["SU3Eps"] = config_dict[""]
            json_dict["flowEpsilon"] = config_dict["flowEpsilon"]
            json_dict["metropolisSeed"] = config_dict["metropolisSeed"]
            json_dict["randomMatrixSeed"] = config_dict["randomMatrixSeed"]
            json.dump(json.JSONEncoder(json_dict),json_file)

    def submitJob(self, job_configurations, system, partition,excluded_nodes=False):
        if excluded_nodes:
            sbatch_exclusions = "#SBATCH --exclude=%s" % excluded_nodes
        else:
            sbatch_exclusions = ""

        for job_config in job_configurations:
            # Retrieving config contents
            binary_filename     = job_config["bin_fn"]
            runName             = job_config["runName"]
            threads             = job_config["threads"]
            beta                = job_config["beta"]
            NSpatial            = job_config["N"]
            NTemporal           = job_config["NT"]
            NTherm              = job_config["NTherm"]
            NCor                = job_config["NCor"] 
            NCf                 = job_config["NCf"]
            NFlows              = job_config["NFlows"]
            NUpdates            = job_config["NUpdates"]
            SU3Eps              = job_config["SU3Eps"]
            flowEpsilon         = job_config["flowEpsilon"]
            storeCfgs           = job_config["storeCfgs"]
            storeThermCfgs      = job_config["storeThermCfgs"]
            hotStart            = job_config["hotStart"]
            subDims             = job_config["subDims"]
            verboseRun          = job_config["verboseRun"]
            uTest               = job_config["uTest"]
            uTestVerbose        = job_config["uTestVerbose"]
            cpu_approx_runtime_hr  = job_config["cpu_approx_runtime_hr"]
            cpu_approx_runtime_min = job_config["cpu_approx_runtime_min"]
            
            # Checks that binary file exists in expected location
            if not os.path.isfile("%s/%s" % (self.CURRENT_PATH,binary_filename)):
                exit("Error: binary file path not in expected location %s/%s" % (self.CURRENT_PATH,binary_filename))

            checkSubDimViability(subDims)

            # Setting job name before creating content file.
            job_name = "{0:<3.2f}beta_{1:<d}cube{2:<d}_{3:<d}threads".format(beta,NSpatial,NTemporal,threads)

            # Setting approximated run time
            estimated_time = "{0:0>2d}:{1:0>2d}:00".format(cpu_approx_runtime_hr,cpu_approx_runtime_min)

            # Setting run-command
            run_command = "mpirun -n {0:<s} {1:<s} {2:<s}".format(threads,binary_filename,json_file)

            # Chosing system
            if system == "smaug":
                # Smaug batch file.
            # #SBATCH --exclude=smaug-b2 excludes an unstable node
                content ='''#!/bin/bash
#SBATCH --job-name={0:<s}
#SBATCH --partition={1:<s}
#SBATCH --ntasks={2:<d}
#SBATCH --time={3:<s}
{4:<s}
{5:<s}
'''.format(job_name, partition, threads, estimated_time, sbatch_exclusions, run_command)

            elif system == "abel":
                # Abel specific commands
                cpu_memory = 3800
                account_name = "nn2977k"
                tasks_per_node = 16 # Maximum number of threads per node
                if threads > tasks_per_node:
                    nodes = threads / tasks_per_node
                    if threads % tasks_per_node != 0:
                        raise ValueError("Tasks(number of threads) have to be divisible by 16.")

                content ='''#!/bin/bash
#SBATCH --job-name={0:<s}
#SBATCH --account={1:<s}
#SBATCH --time={2:<s}
#SBATCH --mem-per-cpu={3:<4d}
#SBATCH --partition={4:<s}
#SBATCH --ntasks={5:<d}
#SBATCH --nodes={6:<1d}
#SBATCH --ntasks-per-node={7:<1d}
{8:<s}

source /cluster/bin/jobsetup

module purge                # clear any inherited modules
module load openmpi.gnu     # loads mpi

set -o errexit               # exit on errors

{9:<s}

'''.format(job_name, account_name, estimated_time, cpu_memory, partition, threads, nodes, tasks_per_node, sbatch_exclusions, run_command)

            job = 'jobfile.slurm'

            # Writies slurm file to be submitted
            if self.dryrun:
                print "Writing %s to slurm batch file: \n\n" % job, content, "\n"
            else:
                outfile = open(job, 'w')
                outfile.write(content)
                outfile.close()

            cmd = ['sbatch', os.getcwd() + "/" + job]

            # Submits job
            if self.dryrun:
                print "> %s %s" % tuple(cmd)
                ID = 0
            else:
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                tmp = proc.stdout.read()
                try:
                    ID = int(tmp.split()[-1])               # ID of job
                except IndexError:
                    print "ERROR: IndexError for line: \n", tmp, "--> exiting", exit(0)

            # Stores job in job dictionary
            self.jobs[ID] = [partition,runName,beta,NSpatial,NTemporal,NCf,NTherm,NCor,NUpdates,NFlows,SU3Eps,threads,bool(storeCfgs),bool(storeThermCfgs),bool(hotStart),' '.join(map(str,subDims)),cpu_approx_runtime_hr,cpu_approx_runtime_min]
            
            # Changes name of job script
            if self.dryrun:
                print '> mv %s job_%d.sh' % (job, ID)
            else:
                os.system('mv %s job_%d.sh' % (job, ID))

        # Updates ID file
        self.updateIdFile()

    def updateIdFile(self):
        # Rewrites the .ids.txt file with additional job ID
        if self.dryrun:
            print "Updating %s file." % self.idFilesName
        else:
            with open(self.idFilesName,"w") as f:
                json.dump(json.JSONEncoder(self.jobs),f)

    def cancelJob(self, jobID):
        if self.dryrun:
            print "scancel %d" % jobID
        else:
            os.system("scancel %d" % jobID)

    def cancelAllJobs(self):
        if self.dryrun:
            for i in self.jobs:
                print "scancel %d" % i
        else:
            for i in self.jobs:
                os.system("scancel %d" % i)

    def showIDwithNb(self):
        header_labels = ['ID', 'Partition', 'Run-name', 'beta', 'N', 'NT', 'NCf', 'NTherm', 'NCor', 'NUpdates', 'NFlows', 'SU3Eps', 'threads', 'storeCfgs', 'storeThermCfgs', 'hotStart', 'subDims', 'Ap.Time[hr]', 'Ap.Time[min]']
        widthChoser = lambda s: 7 if len(s) <= 6 else len(s)+4
        colWidths = [widthChoser(i) for i in header_labels]
        colWidths[2] = 20
        for label,colWidth in zip(header_labels,colWidths):
            print '{0:<{w}}'.format(label, w=colWidth),
        print
        for i in sorted(self.jobs):
            print '{0:<{w}}'.format(i, w=colWidths[0]),
            for j in xrange(len(self.jobs[i])):
                print '{0:<{w}}'.format(self.jobs[i][j], w=colWidths[j+1]),
            print

    def clearIdFile(self):
        self.jobs = {}
        self.updateIdFile()

#------------------------------------------------------------------------------#

def main(args):
    # Default configuration file.
    if not os.path.isdir(os.getcwd() + "/build"):
        raise EnvironmentError("Build folder is not present at location %s." % (os.getcwd() + "/build"))

    # Default config
    config_default = {  "bin_fn"        : "%s/build/GluonicLQCD" % os.getcwd(),
                        "runName"       : "defaultTestRun",
                        "beta"          : 6.0,
                        "N"             : 24,
                        "NT"            : 48,
                        "NTherm"        : 200,
                        "NCor"          : 20,
                        "NCf"           : 100,
                        "NFlows"        : 0,
                        "NUpdates"      : 10,
                        "SU3Eps"        : 0.24,
                        "flowEpsilon"   : 0.01,
                        "threads"       : 64,
                        "storeCfgs"     : 1,
                        "storeThermCfgs": 0,
                        "hotStart"      : 0,
                        "subDims"       : [],
                        "verboseRun"    : 1,
                        "uTest"         : 0,
                        "uTestVerbose"  : 0,
                        "cpu_approx_runtime_hr": 2,
                        "cpu_approx_runtime_min": 0,
                        "cpu_memory"    : 3800,
                        "account_name"  : "nn2977k"}

    ######## Initiating command line parser ########
    description_string = '''
    Program for starting large parallel Lattice Quantum Chromo Dynamics jobs.
    '''
    parser = argparse.ArgumentParser(prog='GluonicLQCD job creator', description=description_string)

    ######## Prints program version if prompted ########
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.1')
    parser.add_argument('--dryrun', default=False, action='store_true', help='Dryrun to no perform any critical actions.')

    ######## Sets up subparsers ########
    subparser = parser.add_subparsers(dest='subparser')

    ######## Job control ########
    sbatch_parser = subparser.add_parser('sbatch', help='Views, stops, clears and list jobs.')
    sbatch_group = sbatch_parser.add_mutually_exclusive_group(required=True)
    sbatch_group.add_argument('--scancel',                      default=False,      type=int,help='Cancel a job of given ID.')
    sbatch_group.add_argument('--scancel_all',                  default=False,      action='store_true',help='Cancel all jobs')
    sbatch_group.add_argument('--list_jobs',                    default=False,      action='store_true',help='List all jobs currently running.')
    sbatch_group.add_argument('--clearIDFile',                  default=False,      action='store_true',help='Clears the job ID file.')

    ######## Manual job setup ########
    job_parser = subparser.add_parser('setup', help='Sets up the job.')
    job_parser.add_argument('system',                           default=False,      type=str, choices=['smaug','abel'],help='Specify system we are running on.')
    job_parser.add_argument('-rn',  '--run_name',               default='run',      type=str,help='Specifiy the run name')
    job_parser.add_argument('-p',   '--partition',              default="normal",   type=str,help='Specify partition to run program on.')
    job_parser.add_argument('-t',   '--threads',                default=False,      type=int,help='Number of threads to run on')
    job_parser.add_argument('-N',   '--NSpatial',               default=False,      type=int,help='spatial lattice dimension')
    job_parser.add_argument('-NT',  '--NTemporal',              default=False,      type=int,help='temporal lattice dimension')
    job_parser.add_argument('-NTh', '--NTherm',                 default=-1,         type=int,help='number of thermalization steps')
    job_parser.add_argument('-NUp', '--NUpdates',               default=False,      type=int,help='number of updates per link')
    job_parser.add_argument('-NCf', '--NConfigs',               default=False,      type=int,help='number of configurations to generate')
    job_parser.add_argument('-NCor', '--NCor',                  default=False,      type=int,help='number of correlation updates to perform')
    job_parser.add_argument('-NFLows','--NFlows',               default=0,          type=int,help='number of flows to perform per configuration')
    job_parser.add_argument('-b',   '--beta',                   default=False,      type=float,help='beta value')
    job_parser.add_argument('-SU3', '--SU3Eps',                 default=False,      type=float,help='SU3 value')
    job_parser.add_argument('-hs', '--hotStart',                default=False,      type=bool,help='Hot start or cold start')
    job_parser.add_argument('-sd', '--subDims',                 default=False,      type=int,nargs=4,help='List of sub lattice dimension sizes, length 4')
    job_parser.add_argument('-sq', '--square',                  default=False,      action='store_true',help='Enforce square sub lattices(or as close as possible).')
    job_parser.add_argument('-sc','--storeCfgs',                default=True,       type=bool,help='Specifying if we are to store configurations')
    job_parser.add_argument('-st', '--storeThermCfgs',          default=False,      type=bool,help='Specifies if we are to store the thermalization plaquettes')
    job_parser.add_argument('-chr', '--cpu_approx_runtime_hr',  default=-1,         type=int,help='Approximate cpu time in hours that will be used')
    job_parser.add_argument('-cmin', '--cpu_approx_runtime_min',default=-1,         type=int,help='Approximate cpu time in minutes that will be used')
    job_parser.add_argument('-ex','--exclude',                  default=False,      type=str,nargs='+',help='Nodes to exclude.')
    
    ######## Abel specific commands ########
    job_parser.add_argument('--cpu_memory',                     default=False,      type=int,help='CPU memory to be allocated to each core')
    job_parser.add_argument('--account_name',                   default=False,      type=str,help='Account name associated to the abel cluster')

    ######## Job load parser ########
    load_parser = subparser.add_parser('load', help='Loads a configuration file into the program')
    load_parser.add_argument('file',                            default=False,      type=str, nargs='+', help='Loads config file')
    load_parser.add_argument('-s','--system',                   default=False,      type=str, required=True,choices=['smaug','abel'],help='Cluster name')
    load_parser.add_argument('-p','--partition',                default="normal",   type=str, help='Partition to run on')

    ######## Unit test parser ########
    unit_test_parser = subparser.add_parser('utest', help='Runs unit tests embedded in the GluonicLQCD program. Will exit when complete.')
    unit_test_parser.add_argument('system',                     default=False,      type=str, choices=['smaug','abel'],help='Specify system we are running on.')
    unit_test_parser.add_argument('-v', '--verbose',            default=False,      action='store_true', help='Prints more information during testing.')

    args = parser.parse_args()
    # args = parser.parse_args(['python', 'makeJobs.py', 'load', 'config_folder/size_scaling_configs/config_16cube32.py', 'config_folder/size_scaling_configs/config_24cube48.py', 'config_folder/size_scaling_configs/config_28cube56.py', 'config_folder/size_scaling_configs/config_32cube64.py', '-s', 'abel'])
    # args = parser.parse_args(["--dryrun","setup","smaug","-ex","smaug-a[1-8]","smaug-b[1-8]","-sq"])
    # args = parser.parse_args(["--dryrun","load","config_folder/test_config_file.py","abel"])
    # args = parser.parse_args(["--dryrun","setup","abel","-rn","test","-subN","4","4","4","4"])

    # Retrieves dryrun bool
    dryrun = args.dryrun

    # Initiates Slurm class for running jobs ect
    s = Slurm(dryrun)

    # Loads a configuration into the world(or multiple!)
    if args.subparser == 'load':
        configurations = [eval(open(load_argument,"r").read()) for load_argument in args.file]
        s.submitJob(configurations,args.system,args.partition)
    elif args.subparser == 'setup':
        # Note: with setup, can only submit a single job at the time
        partition = "normal"
        excluded_nodes = ""
        system = args.system
        if args.run_name:
            config_default["runName"] = args.run_name
        if args.partition:
            partition = args.partition
        if args.threads:
            config_default["threads"] = args.threads
        if args.NSpatial:
            config_default["N"] = args.NSpatial
        if args.NTemporal:
            config_default["NT"] = args.NTemporal
        if args.NTherm != -1:
            config_default["NTherm"] = args.NTherm
        if args.NUpdates:
            config_default["NUpdates"] = args.NUpdates
        if args.NCor:
            config_default["NCor"] = args.NCor
        if args.NConfigs:
            config_default["NCf"] = args.NConfigs
        if args.NFlows:
            config_default["NFlows"] = args.NFlows
        if args.beta:
            config_default["beta"] = args.beta
        if args.SU3Eps:
            config_default["SU3Eps"] = args.SU3Eps
        if args.hotStart:
            config_default["hotStart"] = int(args.hotStart)
        if args.subDims:
            checkSubDimViability(args.subDims)
            config_default["subDims"] = args.subDims
        if args.square:
            config_default["subDims"] = createSquare(config_default["threads"],config_default["N"],config_default["NT"])
        if args.storeCfgs:
            config_default["storeCfgs"] = args.storeCfgs
        if args.exclude:
            excluded_nodes = ','.join(args.exclude)
        if args.cpu_approx_runtime_hr != -1 and system == "abel":
            config_default["cpu_approx_runtime_hr"] = args.cpu_approx_runtime_hr
        if args.cpu_approx_runtime_min != -1 and system == "abel":
            config_default["cpu_approx_runtime_min"] = args.cpu_approx_runtime_min
        if args.account_name and system == "abel":
            config_default["account_name"] = args.account_name
        # Submitting job
        s.submitJob([config_default],system,partition,excluded_nodes)
    elif args.subparser == 'sbatch':
        if args.scancel:
            s.cancelJob(args.scancel)
        if args.scancel_all:
            s.cancelAllJobs()
        if args.list_jobs:
            s.showIDwithNb()
        if args.clearIDFile:
            s.clearIdFile()
    elif args.subparser == 'utest':
        config_default["uTest"] = 1
        config_default["cpu_approx_runtime_hr"] = 0
        config_default["cpu_approx_runtime_min"] = 10
        partition = "normal"
        excluded_nodes = ""
        system = args.system
        if args.verbose:
            config_default["uTestVerbose"] = int(args.verbose)
        s.submitJob([config_default],system,partition,excluded_nodes)
    else:
        print 'Parse error: %s \n--> exiting' % args
        exit(0)

if __name__ == '__main__':
    main(sys.argv[1:])