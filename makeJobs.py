import os, subprocess, time, sys

'''
TODO: 
[ ] Creat command-line tools!
[ ] Add Abel usability and Smaug usability.
[ ] make command line such that I can show .ids.txt
[ ] Add pwd in path to avoid any confusion
[ ] Create folders needed
'''

class Slurm:
    def __init__(self):
        self.idFilesName = '.ids.txt'
        if os.path.isfile(self.idFilesName):
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
            for dim in subDims:
                if dim <= 2: exit("Error: %d is not a valid dimension" % dim)

            content ='''
#!/bin/bash
#SBATCH --job-name=%s_Nb%d
#SBATCH --account=nn2977k
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=3800M
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
source /cluster/bin/jobsetup
module purge
module load openmpi.gnu
set -o errexit
cp build/GluonicLQCD $SCRATCH
chkfile output.txt
cd $SCRATCH
./CCD "%s" %d %d %.02f %.02e %d %d
            '''

            # #SBATCH --exclude=smaug-b2 excludes an unstable node
            content ='''#!/bin/bash
#SBATCH --partition=%s
#SBATCH --ntasks=%d
#SBATCH --time=12:00:00
#SBATCH --job-name=%3.2fbeta_%dcube%d_%dthreads
mpirun -n %d %s %s %d %d %d %d %d %d %.2f %.2f %1d %1d %1d %s
        '''%(partition,threads,beta,NSpatial,NTemporal,threads,threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps,storeCfgs,storeThermCfgs,hotStart,' '.join(map(str,subDims)))
            job = 'jobfile.slurm'
            outfile  = open(job, 'w')
            outfile.write(content)
            outfile.close()
            cmd = ['sbatch', os.getcwd() + "/" + job] # Do I need the os.getcwd()?
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
            ID = int(tmp.split()[-1])               # ID of job
            self.jobs[ID] = [partition,runName,beta,NSpatial,NTemporal,NCf,NTherm,NCor,NUpdates,SU3Eps,threads,bool(storeCfgs),bool(storeThermCfgs),bool(hotStart),' '.join(map(str,subDims))]
            os.system('mv %s job_%d.sh'%(job, ID))  # Change name of jobScript
        self.updateIdFile()

    def updateIdFile(self):
        f = open(self.idFilesName,"w")
        f.write(str(self.jobs))
        f.close()

    def cancelJob(self, jobID):
        os.system('scancel %d'%jobID)

    def cancelAllJobs(self):
        for i in self.jobs:
            os.system('scancel %d'%i)

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

if __name__ == '__main__':
    # Variables
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

    s = Slurm()
    s.submitJob([job_configuration1], partition)
    #s.cancelAllJobs()
    s.showIDwithNb()
