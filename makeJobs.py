import os, subprocess, time, sys, argparse

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
    return N

def checkSubDimViability(subDims):
    if len(subDims) != 4 and sum([type(i) == int for i in subDims]) != 4:
        raise ValueError("%g is not a valid set of sub dimensions." % subDims)

class Slurm:
    def __init__(self, dryrun):
        self.dryrun = dryrun
        self.CURRENT_PATH = os.getcwd()
        # Checking that we have an output folder.
        if not os.path.isdir('%s/output' % self.CURRENT_PATH):
            if self.dryrun:
                print '> mkdir %s/output' % self.CURRENT_PATH
            else:
                os.mkdir('%s/output' % self.CURRENT_PATH)

        # Checking the .ids.txt file for previous jobs and storing them in job-list.
        self.idFilesName = '.ids.txt'
        if os.path.isfile(self.idFilesName) and not self.dryrun:
            self.jobs =  eval(open(self.idFilesName,"r").read())
        else:
            self.jobs = {}

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
            NUpdates            = job_config["NUpdates"]
            SU3Eps              = job_config["SU3Eps"]
            storeCfgs           = job_config["storeCfgs"]
            storeThermCfgs      = job_config["storeThermCfgs"]
            hotStart            = job_config["hotStart"]
            subDims             = job_config["subDims"]
            cpu_approx_runtime  = job_config["cpu_approx_runtime"]
            
            # Error catching before submitting job is nice. MIGHT ADD MORE
            for dim in subDims:
                if dim <= 2: exit("Error: %d is not a valid dimension" % dim)

            # Chosing system
            if system == "smaug":
                # Smaug batch file.
            # #SBATCH --exclude=smaug-b2 excludes an unstable node
                content ='''#!/bin/bash
#SBATCH --job-name={2:<3.2f}beta_{3:<d}cube{4:<d}_{5:<d}threads
#SBATCH --partition={0:<s}
#SBATCH --ntasks={1:<d}
#SBATCH --time={21:0>2d}:00:00
{23:<s}
mpirun -n {6:<d} {7:<s} {8:<s} {9:<d} {10:<d} {11:<d} {12:<d} {13:<d} {14:<d} {15:<.2f} {16:<.2f} {17:<1d} {18:<1d} {19:<1d} {22:<s} {20:<s}
'''.format( partition,threads,beta,NSpatial,NTemporal,threads,
                threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)),
                cpu_approx_runtime,self.CURRENT_PATH,sbatch_exclusions)
            elif system == "abel":
                # cpu_memory = job_config["cpu_memory"]
                # account_name = job_config["account_name"] # Make system specific?
                # nodes = job_config["nodes"]
                # tasks_per_node = job_config["tasks_per_node"]

                # Abel specific commands
                cpu_memory = 3800
                account_name = "nn2977k"
                tasks_per_node = 16 # Maximum number of threads per node
                if threads > tasks_per_node:
                    nodes = threads / tasks_per_node
                    if threads % tasks_per_node != 0:
                        raise ValueError("Tasks(number of threads) have to be divisible by 16.")

                # Bash file to run on Abel
                content ='''#!/bin/bash
#SBATCH --job-name={0:<3.2f}beta_{1:<d}cube{2:<d}_{3:<d}threads
#SBATCH --account={23:<s}
#SBATCH --time={21:0>2d}:00:00
#SBATCH --mem-per-cpu={22:<4d}M
#SBATCH --partition={4:<s}
#SBATCH --ntasks={5:<d}
#SBATCH --nodes={24:<1d}
#SBATCH --ntasks-per-node={25:<1d}
{27:<s}

source /cluster/bin/jobsetup

module purge                # clear any inherited modules
module load openmpi.gnu     # loads mpi

set -o errexit               # exit on errors

#chkfile output
#chkfile "*.bin"
#chkfile "*.dat"
#chkfile output.txt

#cp build/GluonicLQCD $SCRATCH
#cd $SCRATCH
#mkdir output

mpirun -n {6:<d} {7:<s} {8:<s} {9:<d} {10:<d} {11:<d} {12:<d} {13:<d} {14:<d} {15:<.2f} {16:<.2f} {17:<1d} {18:<1d} {19:<1d} {26:<s} {20:<s}
'''.format(beta,NSpatial,NTemporal,threads,partition,threads,
                threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)),
                cpu_approx_runtime,cpu_memory,account_name,nodes,tasks_per_node,self.CURRENT_PATH,sbatch_exclusions)

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
            self.jobs[ID] = [partition,runName,beta,NSpatial,NTemporal,NCf,NTherm,NCor,NUpdates,SU3Eps,threads,bool(storeCfgs),bool(storeThermCfgs),bool(hotStart),' '.join(map(str,subDims)),cpu_approx_runtime]
            
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
            f = open(self.idFilesName,"w")
            f.write(str(self.jobs))
            f.close()

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
        header_labels = ['ID', 'Partition', 'Run-name', 'beta', 'N', 'NT', 'NCf', 'NTherm', 'NCor', 'NUpdates', 'SU3Eps', 'threads', 'storeCfgs', 'storeThermCfgs', 'hotStart', 'subDims', 'CPU-estimate[hours]']
        widthChoser = lambda s: 7 if len(s) <= 6 else len(s)+4
        colWidths = [widthChoser(i) for i in header_labels]
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
                        "NUpdates"      : 10,
                        "SU3Eps"        : 0.24,
                        "threads"       : 64,
                        "storeCfgs"     : 1,
                        "storeThermCfgs": 0,
                        "hotStart"      : 0,
                        "subDims"       : [],
                        "cpu_approx_runtime": 2,
                        "cpu_memory"    : 3800,
                        "account_name"  : "nn2977k"}

    # Initiating command line parser.
    description_string = '''
    Program for starting large parallel Lattice Quantum Chromo Dynamics jobs.
    '''
    parser = argparse.ArgumentParser(prog='GluonicLQCD job creator', description=description_string)

    # Prints program version if prompted
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.1')
    parser.add_argument('--dryrun', default=False, action='store_true', help='Dryrun to no perform any critical actions.')

    subparser = parser.add_subparsers(dest='subparser')

    # Job control
    sbatch_parser = subparser.add_parser('sbatch', help='Views, stops, clears and list jobs.')
    sbatch_group = sbatch_parser.add_mutually_exclusive_group(required=True)
    sbatch_group.add_argument('--scancel',          default=False,      type=int,help='Cancel a job of given ID.')
    sbatch_group.add_argument('--scancel_all',      default=False,      action='store_true',help='Cancel all jobs')
    sbatch_group.add_argument('--list_jobs',        default=False,      action='store_true',help='List all jobs currently running.')
    sbatch_group.add_argument('--clearIDFile',      default=False,      action='store_true',help='Clears the job ID file.')

    # Job setup
    job_parser = subparser.add_parser('setup', help='Sets up the job.')
    job_parser.add_argument('system',               default=False,      type=str, choices=['smaug','abel'],help='Specify system we are running on.')
    job_parser.add_argument('-rn',  '--run_name',   default='run',      type=str,help='Specifiy the run name')
    job_parser.add_argument('-p',   '--partition',  default="normal",   type=str,help='Specify partition to run program on.')
    job_parser.add_argument('-t',   '--threads',    default=False,      type=int,help='Number of threads to run on')
    job_parser.add_argument('-N',   '--NSpatial',   default=False,      type=int,help='spatial lattice dimension')
    job_parser.add_argument('-NT',  '--NTemporal',  default=False,      type=int,help='temporal lattice dimension')
    job_parser.add_argument('-NTh', '--NTherm',     default=False,      type=int,help='number of thermalization steps')
    job_parser.add_argument('-NUp', '--NUpdates',   default=False,      type=int,help='number of updates per link')
    job_parser.add_argument('-NCf', '--NConfigs',   default=False,      type=int,help='number of configurations to generate')
    job_parser.add_argument('-b',   '--beta',       default=False,      type=float,help='beta value')
    job_parser.add_argument('-SU3', '--SU3Eps',     default=False,      type=float,help='SU3 value')
    job_parser.add_argument('-hs', '--hotStart',    default=False,      type=bool,help='Hot start or cold start')
    job_parser.add_argument('-sd', '--subDims',     default=False,      type=int,nargs=4,help='List of sub lattice dimension sizes, length 4')
    job_parser.add_argument('-sq', '--square',      default=False,      action='store_true',help='Enforce square sub lattices(or as close as possible).')
    job_parser.add_argument('-sc','--storeCfgs',    default=True,       type=bool,help='Specifying if we are to store configurations')
    job_parser.add_argument('-st', '--storeThermCfgs',    default=False,type=bool,help='Specifies if we are to store the thermalization plaquettes')
    job_parser.add_argument('-c', '--cpu_approx_runtime', default=False,type=int,help='Approximate cpu time that will be used')
    job_parser.add_argument('-ex','--exclude',      default=False,      type=str,nargs='+',help='Nodes to exclude.')
    # Only to be used at abel
    job_parser.add_argument('--cpu_memory',         default=False,      type=int,help='CPU memory to be allocated to each core')
    job_parser.add_argument('--account_name',       default=False,      type=str,help='Account name associated to the abel cluster')

    # Job load parser
    load_parser = subparser.add_parser('load', help='Loads a configuration file into the program')
    load_parser.add_argument('file',                default=False,      type=str, nargs='+', help='Loads config file')
    load_parser.add_argument('-s','--system',       default=False,      type=str, required=True,choices=['smaug','abel'],help='Cluster name')
    load_parser.add_argument('-p','--partition',    default="normal",   type=str, help='Partition to run on')

    args = parser.parse_args()
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
        if args.NTherm:
            config_default["NTherm"] = args.NTherm
        if args.NUpdates:
            config_default["NUpdates"] = args.NUpdates
        if args.NConfigs:
            config_default["NCf"] = args.NConfigs
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
        if args.cpu_approx_runtime and system == "abel":
            config_default["cpu_approx_runtime"] = args.cpu_approx_runtime
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
    else:
        print 'Parse error: %s \n--> exiting' % args
        exit(0)

if __name__ == '__main__':
    """
    TODO:
    """
    main(sys.argv[1:])