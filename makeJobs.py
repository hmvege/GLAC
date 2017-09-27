import os, subprocess, time

class Slurm:
    def __init__(self, partition):
        self.partition = partition
        self.idFilesName = '.ids.txt'
        if os.path.isfile(self.idFilesName):
            self.jobs =  eval(open(self.idFilesName,"r").read())
        else:
            self.jobs = {}

    def submitJob(self, job_configurations):
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

            content ='''#!/bin/bash
#SBATCH --partition=%s
#SBATCH --ntasks=64  
#SBATCH --time=01:00:00
#SBATCH --job-name=%3.2fbeta_%dcube%d_%dthreads
mpirun -n %d %s %s %d %d %d %d %d %d %.2f %.2f
        '''%(self.partition,beta,NSpatial,NTemporal,threads,threads,binary_filename,runName,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps)
            job = 'jobfile.slurm'
            outfile  = open(job, 'w')
            outfile.write(content)
            outfile.close()
            cmd = ['sbatch', os.getcwd() + "/" + job] # Do I need the os.getcwd()?
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
            ID = int(tmp.split()[-1])               # ID of job
            self.jobs[ID] = [runName,beta,NSpatial,NTemporal,NCf,NTherm,NCor,NUpdates,SU3Eps,threads,threads]
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
        colWidth = 12
        print '{0:<{w}} {1:<{w}} {2:<{w}} {3:<{w}} {4:<{w}} {5:<{w}} {6:<{w}} {7:<{w}}'.format('ID', 'Run-name', 'beta', 'N', 'NT', 'NCf', 'NTherm', 'NCor', 'NUpdates', 'SU3Eps', 'threads', w=colWidth)
        for i in sorted(self.jobs):
            print '{0:<{w}}'.format(i, w=colWidth),
            for j in xrange(len(self.jobs[i])):
                print '{0:<{w}}'.format(self.jobs[i][j], w=colWidth),
            print

    def clearIdFile(self):
        self.jobs = {}
        self.updateIdFile()

#------------------------------------------------------------------------------#

if __name__ == '__main__':
    # Variables
    job_configuration1 = {  "bin_fn"    : "%s/build/GluonicLQCD" % os.getcwd(),
                            "runName"   : "testRun1",
                            "beta"      : 6.0,
                            "N"         : 24,
                            "NT"        : 48,
                            "NTherm"    : 200,
                            "NCor"      : 20,
                            "NCf"       : 100,
                            "NUpdates"  : 10,
                            "SU3Eps"    : 0.24,
                            "threads"   : 64}

    job_configuration2 = {}

    s = Slurm("normal")
    s.submitJob([job_configuration1])
    #s.cancelAllJobs()
    s.showIDwithNb()
