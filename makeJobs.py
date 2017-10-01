import os, subprocess, time, sys, argparse

'''
TODO: 
[ ] Create command-line tools!
[ ] Add Abel usability and Smaug usability.
[ ] make command line such that I can show .ids.txt
[x] Add pwd in path to avoid any confusion
[ ] Create folders needed
[ ] Fix specifying dimensions (C++).
[ ] Fix writing out to file functions (C++).
'''

class Slurm:
    def __init__(self, dryrun):
        self.dryrun = dryrun

        # Checking that we have an output folder.
        if not os.path.isdir('%s/output' % os.getcwd()):
            if self.dryrun:
                print '> mkdir %s/output' % os.getcwd()
            else:
                os.mkdir('%s/output' % os.getcwd())

        # Checking the .ids.txt file for previous jobs and storing them in job-list.
        self.idFilesName = '.ids.txt'
        if os.path.isfile(self.idFilesName) and not self.dryrun:
            self.jobs =  eval(open(self.idFilesName,"r").read())
        else:
            self.jobs = {}

    def submitJob(self, system, job_configurations, partition):
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
            
            # Error catching before submitting job is nice.
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
mpirun -n {6:<d} {7:<s} {8:<s} {9:<d} {10:<d} {11:<d} {12:<d} {13:<d} {14:<d} {15:<.2f} {16:<.2f} {17:<1d} {18:<1d} {19:<1d} {20:<s}
'''.format( partition,threads,beta,NSpatial,NTemporal,threads,
                threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)),
                cpu_approx_runtime)
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
                    nodes = tasks / tasks_per_node
                    if tasks % tasks_per_node != 0:
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

mpirun -n {6:<d} {7:<s} {8:<s} {9:<d} {10:<d} {11:<d} {12:<d} {13:<d} {14:<d} {15:<.2f} {16:<.2f} {17:<1d} {18:<1d} {19:<1d} {20:<s}
'''.format(beta,NSpatial,NTemporal,threads,partition,threads,
                threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)),
                cpu_approx_runtime,cpu_memory,account_name,nodes,tasks_per_node)

            job = 'jobfile.slurm'

            # Writies slurm file to be submitted
            if self.dryrun:
                print "Writing %s to slurm batch file: \n\n" % job, content, "\n"
            else:
                outfile  = open(job, 'w')
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
    description_string = '''
    Program for starting large parallel Lattice Quantum Chromo Dynamics jobs.

    TODO: 
    '''
    parser = argparse.ArgumentParser(prog='GluonicLQCD job creator', description=description_string)

    # Prints program version if prompted
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.1')
    parser.add_argument('--dryrun', default=False, action='store_true',help='Dryrun to no perform any critical actions.')

    subparser = parser.add_subparsers(dest='subparser')

    # Job control
    sbatch_parser = subparser.add_parser('sbatch', help='Views, stops, clears and list jobs.')
    sbatch_group = sbatch_parser.add_mutually_exclusive_group(required=True)
    sbatch_group.add_argument('--scancel',    default=False,      type=int,help='Cancel a job of given ID.')
    sbatch_group.add_argument('--scancel_all',default=False,      action='store_true',help='Cancel all jobs')
    sbatch_group.add_argument('--list_jobs',  default=False,      action='store_true',help='List all jobs currently running.')
    sbatch_group.add_argument('--clearIDFile',default=False,      action='store_true',help='Clears the job ID file.')

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
    job_parser.add_argument('--SU3',                default=False,      type=float,help='SU3 value')

    args = parser.parse_args()

    # Retrieves dryrun bool
    dryrun = args.dryrun

    # Initiates Slurm class for running jobs ect
    s = Slurm(dryrun)

    print args
    if args.subparser == 'setup':
        s.submitJob()
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
        print 'Parse error:', args
        exit(0)

if __name__ == '__main__':
    main(sys.argv[1:])
    # Variables
    # exit(0)
    # system = "abel"
    # dryrun = True

    # if len(sys.argv) > 1:
    #     partition = str(sys.argv[1])
    # else:
    #     partition = "normal"

    # # sub_lattice_dimensions = [12,12,12,6]
    # default_configuration = {  "bin_fn"        : "%s/build/GluonicLQCD" % os.getcwd(),
    #                         "runName"       : "testRun1",
    #                         "beta"          : 6.0,
    #                         "N"             : 24,
    #                         "NT"            : 48,
    #                         "NTherm"        : 200,
    #                         "NCor"          : 20,
    #                         "NCf"           : 100,
    #                         "NUpdates"      : 10,
    #                         "SU3Eps"        : 0.24,
    #                         "threads"       : 64,
    #                         "storeCfgs"     : 1,
    #                         "storeThermCfgs": 0,
    #                         "hotStart"      : 0,
    #                         "subDims"       : [12,12,12,6],
    #                         "cpu_approx_runtime": 2,
    #                         "cpu_memory"    : 3800,
    #                         "account_name"  : "nn2977k",
    #                         "nodes"         : 4,
    #                         "tasks_per_node": 16}
    # job_configuration2 = {}

    # s = Slurm(dryrun)
    # s.submitJob(system, [job_configuration1], partition)
    # # s.cancelAllJobs()
    # # s.clearIdFile()
    # # s.cancelJob(ID)
    # s.showIDwithNb()
