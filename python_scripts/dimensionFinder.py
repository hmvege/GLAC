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
		NSpatial = 24
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
	m_N = [10, 10, 10]
	i = 9
	j = 9
	k = 9
	l = 9
	print i + m_N[0]*(j + m_N[1]*(k + m_N[2]*l))
	main(sys.argv[1:])