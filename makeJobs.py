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
            
            # Error catching before submitting job is nice.
            for dim in subDims:
                if dim <= 2: exit("Error: %d is not a valid dimension" % dim)

            if self.system == "smaug":
                # Smaug batch file.
            # #SBATCH --exclude=smaug-b2 excludes an unstable node
                content ='''#!/bin/bash
#SBATCH --partition=%s
#SBATCH --ntasks=%d
#SBATCH --time=12:00:00
#SBATCH --job-name=%3.2fbeta_%dcube%d_%dthreads
mpirun -n %d %s %s %d %d %d %d %d %d %.2f %.2f %1d %1d %1d %s
'''%(partition,threads,beta,NSpatial,NTemporal,threads,threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)))
            elif self.system == "abel":
                # Abel batch file.
                content ='''#!/bin/bash
#SBATCH --job-name=%3.2fbeta_%dcube%d_%dthreads
#SBATCH --account=hmvege
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=3800M
#SBATCH --partition=%s
#SBATCH --ntasks=%d
## #SBATCH --cpus-per-task=16

source /cluster/bin/jobsetup

module purge                # clear any inherited modules
module load openmpi.gnu     # loads mpi

set-o errexit               # exit on errors

chkfile output/
chkfile "*.bin"
chkfile "*.dat"
#chkfile output.txt

cp build/GluonicLQCD $SCRATCH
cd $SCRATCH
mkdir output

mpirun -n %d %s %s %d %d %d %d %d %d %.2f %.2f %1d %1d %1d %s
'''%(beta,NSpatial,NTemporal,threads,partition,threads,threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)))

            job = 'jobfile.slurm'

            if self.dryrun:
                print "Writing %s to slurm batch file: \n\n" % job, content, "\n"
            else:
                outfile  = open(job, 'w')
                outfile.write(content)
                outfile.close()

            cmd = ['sbatch', os.getcwd() + "/" + job] # Do I need the os.getcwd()?

            if self.dryrun:
                print "> %s %s" % tuple(cmd)
                ID = 0
            else:
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                tmp = proc.stdout.read()
                ID = int(tmp.split()[-1])               # ID of job

            self.jobs[ID] = [partition,runName,beta,NSpatial,NTemporal,NCf,NTherm,NCor,NUpdates,SU3Eps,threads,bool(storeCfgs),bool(storeThermCfgs),bool(hotStart),' '.join(map(str,subDims))]
            
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
        header_labels = ['ID', 'Partition', 'Run-name', 'beta', 'N', 'NT', 'NCf', 'NTherm', 'NCor', 'NUpdates', 'SU3Eps', 'threads', 'storeCfgs', 'storeThermCfgs', 'hotStart','subDims']
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

# def main(args):
#     parser = argparse.ArgumentParser(prog='GluonicLQCD job creator', description='Program starting Lattice Quantum Chromo Dynamics.')

#     # Prints program version if prompted
#     parser.add_argument('--version', action='version', version='%(prog)s 0.9.2')

#     # Main argument, must have this one
#     parser.add_argument('element',                      default=False,  type=str,   nargs=1,    help='takes the type of element. E.g. He')

#     # Possible choices
#     parser.add_argument('-lf',  '--local_file',         default=None,   type=str,               help='takes a .html file for an atomic spectra from nist.org')
#     parser.add_argument('-fn',  '--filename',           default=None,   type=str,               help='output filename')
#     parser.add_argument('-p',   '--parallel',           default=False,  action='store_const',   const=True, help='enables running in parallel')
#     parser.add_argument('-n',   '--num_processors',     default=4,      type=int,               help='number of processors, default=4')
#     parser.add_argument('-ln',  '--length',             default=10,     type=float,             help='length in seconds, default=10')
#     parser.add_argument('-hz',  '--hertz',              default=440,    type=int,               help='frequency, default=440')
#     parser.add_argument('-amp', '--amplitude',          default=0.01,   type=float,             help='amplitude of track, default=0.01')
#     parser.add_argument('-sr',  '--sampling_rate',      default=44100,  type=int,               help='sampling rate, default=44100')
#     parser.add_argument('-cf',  '--convertion_factor',  default=100,    type=float,             help='factor to pitch-shift spectra by, default=100')
#     parser.add_argument('-wlc', '--wavelength_cutoff',  default=2.5e-1, type=float,             help='inverse wavelength to cutoff lower tones, default=2.5e-1.')
#     parser.add_argument('-bt',  '--beat_cutoff',        default=1e-2,   type=float,             help='removes one wavelength if two wavelengths have |lambda-lambda_0| > beat_cutoff, default=1e-2')

#     args = parser.parse_args()
#     element = args.element[0]
#     if not element_search(element):
#         sys.exit('Element %s not found.' % element)

#     Sound = ElementSound(element, args.local_file, args.filename, args.parallel, args.num_processors)
#     Sound.remove_beat(args.beat_cutoff)
#     Sound.create_sound(args.length, args.hertz, args.amplitude, args.sampling_rate, args.convertion_factor, args.wavelength_cutoff)


if __name__ == '__main__':
    # main(sys.args[1:])
    # Variables
    system = "abel"
    dryrun = False

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
                            "subDims"       : [12,12,12,6]}
    job_configuration2 = {}

    s = Slurm(system, dryrun)
    s.submitJob([job_configuration1], partition)
    # s.cancelAllJobs()
    # s.clearIdFile()
    # s.cancelJob(ID)
    s.showIDwithNb()
