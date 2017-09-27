import os, subprocess, time

class Slurm:
    def __init__(self):
        self.idFilesName = '.ids.txt'
        if os.path.isfile(self.idFilesName):
            self.jobs =  eval(open(self.idFilesName,"r").read())
        else:
            self.jobs = {}

    def submitJob(self, partition, beta, NSpatial, NTemporal, threads, binary_filename):
        for Nb in range(minNb, maxNb):
            content ='''#!/bin/bash
#SBATCH --partition=%s
#SBATCH --ntasks=64  
#SBATCH --time=01:00:00
#SBATCH --job-name=%3.2fbeta_%dcube%d_%dthreads
mpirun -n %d %s %d %d %d %d %d %d %.2f %.2f
        '''%(partition,beta,NSpatial,NTemporal,threads,binary_filename,NSpatial,NTemporal,NTherm,NCor,NCf,NUpdates,beta,SU3Eps)
            job = 'jobfile.slurm'
            outfile  = open(job, 'w')
            outfile.write(content)
            outfile.close()
            cmd = ['sbatch', job]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
            ID = int(tmp.split()[-1])               # ID of job
            self.jobs[ID] = [ threads]
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
        print '{0:<{w}} {1:<{w}} {2:<{w}} {3:<{w}} {4:<{w}} {5:<{w}} {6:<{w}} {7:<{w}}'.format('ID', 'Model', 'Nh', 'Nb', 'var', 'precision', 'degree', 'threads', w=colWidth)
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
    job_configuration = {""}
    
    s = Slurm()
    s.submitJob(10,11)
    #s.cancelAllJobs()
    s.showIDwithNb()
