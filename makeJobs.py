import os, subprocess, time, sys, argparse

'''
TODO: 
[ ] Creat command-line tools!
[ ] Add Abel usability and Smaug usability.
[ ] make command line such that I can show .ids.txt
[ ] Add pwd in path to avoid any confusion
[ ] Create folders needed
'''

class Slurm:
    def __init__(self, system, dryrun):
        self.system = system
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

    def submitJob(self, job_configurations, partition):
        for job_config in job_configurations:
            # Retrieving config contents
            binary_filename = job_config["bin_fn"]
            runName = job_config["runName"]
            threads = job_config["threads"]
            beta = job_config["beta"]
            NSpatial = job_config["N"]
            NTemporal = job_config["NT"]
            NTherm = job_config["NTherm"]
            NCor = job_config["NCor"] 
            NCf = job_config["NCf"]
            NUpdates = job_config["NUpdates"]
            SU3Eps = job_config["SU3Eps"]
            storeCfgs = job_config["storeCfgs"]
            storeThermCfgs = job_config["storeThermCfgs"]
            hotStart = job_config["hotStart"]
            subDims = job_config["subDims"]
            cpu_approx_runtime = job_config["cpu_approx_runtime"]
            # Abel specific items
            cpu_memory = job_config["cpu_memory"]
            account_name = job_config["account_name"] # Make system specific?
            nodes = job_config["nodes"]
            tasks_per_node = job_config["tasks_per_node"]
            
            # Error catching before submitting job is nice.
            for dim in subDims:
                if dim <= 2: exit("Error: %d is not a valid dimension" % dim)

            # Chosing system
            if self.system == "smaug":
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
            elif self.system == "abel":
                # Abel batch file.
                # threads, 
                # if threads > 16:

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

            if self.dryrun:
                print "Writing %s to slurm batch file: \n\n" % job, content, "\n"
            else:
                outfile  = open(job, 'w')
                outfile.write(content)
                outfile.close()

            cmd = ['sbatch', os.getcwd() + "/" + job]

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

            self.jobs[ID] = [partition,runName,beta,NSpatial,NTemporal,NCf,NTherm,NCor,NUpdates,SU3Eps,threads,bool(storeCfgs),bool(storeThermCfgs),bool(hotStart),' '.join(map(str,subDims)),cpu_approx_runtime]
            
            if self.dryrun:
                print '> mv %s job_%d.sh' % (job, ID)
            else:
                os.system('mv %s job_%d.sh' % (job, ID))  # Change name of jobScript
        self.updateIdFile()

    def updateIdFile(self):
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
    parser = argparse.ArgumentParser(prog='GluonicLQCD job creator', description='Program starting Lattice Quantum Chromo Dynamics.')

    # Prints program version if prompted
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.1')

    # # Main argument, must have this one
    # parser.add_argument('element',                      default=False,  type=str,   nargs=1,    help='takes the type of element. E.g. He')

    # # Possible choices
    # parser.add_argument('-lf',  '--local_file',         default=None,   type=str,               help='takes a .html file for an atomic spectra from nist.org')
    # parser.add_argument('-fn',  '--filename',           default=None,   type=str,               help='output filename')
    # parser.add_argument('-p',   '--parallel',           default=False,  action='store_const',   const=True, help='enables running in parallel')
    # parser.add_argument('-n',   '--num_processors',     default=4,      type=int,               help='number of processors, default=4')
    # parser.add_argument('-ln',  '--length',             default=10,     type=float,             help='length in seconds, default=10')
    # parser.add_argument('-hz',  '--hertz',              default=440,    type=int,               help='frequency, default=440')
    # parser.add_argument('-amp', '--amplitude',          default=0.01,   type=float,             help='amplitude of track, default=0.01')
    # parser.add_argument('-sr',  '--sampling_rate',      default=44100,  type=int,               help='sampling rate, default=44100')
    # parser.add_argument('-cf',  '--convertion_factor',  default=100,    type=float,             help='factor to pitch-shift spectra by, default=100')
    # parser.add_argument('-wlc', '--wavelength_cutoff',  default=2.5e-1, type=float,             help='inverse wavelength to cutoff lower tones, default=2.5e-1.')
    # parser.add_argument('-bt',  '--beat_cutoff',        default=1e-2,   type=float,             help='removes one wavelength if two wavelengths have |lambda-lambda_0| > beat_cutoff, default=1e-2')

    # args = parser.parse_args()
    # element = args.element[0]
    # if not element_search(element):
    #     sys.exit('Element %s not found.' % element)

    # Sound = ElementSound(element, args.local_file, args.filename, args.parallel, args.num_processors)
    # Sound.remove_beat(args.beat_cutoff)
    # Sound.create_sound(args.length, args.hertz, args.amplitude, args.sampling_rate, args.convertion_factor, args.wavelength_cutoff)


if __name__ == '__main__':
    # main(sys.argv[1:])
    # Variables
    system = "abel"
    dryrun = True

    if len(sys.argv) > 1:
        partition = str(sys.argv[1])
    else:
        partition = "normal"

    # sub_lattice_dimensions = [12,12,12,6]
    job_configuration1 = {  "bin_fn"        : "%s/build/GluonicLQCD" % os.getcwd(),
                            "runName"       : "testRun1",
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
                            "subDims"       : [12,12,12,6],
                            "cpu_approx_runtime": 2,
                            "cpu_memory"    : 3800,
                            "account_name"  : "nn2977k",
                            "nodes"         : 4,
                            "tasks_per_node": 16}
    job_configuration2 = {}

    s = Slurm(system, dryrun)
    s.submitJob([job_configuration1], partition)
    # s.cancelAllJobs()
    # s.clearIdFile()
    # s.cancelJob(ID)
    s.showIDwithNb()
