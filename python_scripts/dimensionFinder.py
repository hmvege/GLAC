import sys, numpy as np

def getArgMaxIndex(N):
	val = N[0]
	index = 0
	for i in xrange(4):
		if N[i] > val:
			val = N[i]
			index = i
	return index

def main(args):
	if len(args) == 0:
		NSpatial = 16
		NTemporal = 2*NSpatial
		numprocs = 32
	else:
		NSpatial = int(args[0])
		NTemporal = int(args[1])
		numprocs = int(args[2])

	restProc = numprocs;
	N = np.zeros(4)

	for i in xrange(3):
		N[i] = NSpatial
		N[3] = NTemporal;

	while restProc >= 2:
		max_index = getArgMaxIndex(N)
		N[max_index] /= 2
		restProc /= 2
		if (restProc < 2):
			break


	print "PROCESSORS: %d" % numprocs
	print "N:          %d" % NSpatial
	print "NT:         %d" % NTemporal
	print "SUB-LATTICE-DIMENSIONS: ", N[::-1]

if __name__ == '__main__':
	main(sys.argv[1:])