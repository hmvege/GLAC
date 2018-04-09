import numpy as np

def getArgMaxIndex(N):
	val = N[0]
	index = 0
	for i in xrange(4):
		if N[i] > val:
			val = N[i]
			index = i
	return index

def get_sub_dims(N, NT, numprocs):
	restProc = numprocs;
	NSub = np.zeros(4)

	for i in xrange(3):
		NSub[i] = N
	NSub[3] = NT;

	while restProc >= 2:
		max_index = getArgMaxIndex(NSub)
		NSub[max_index] /= 2.
		restProc /= 2.
		if NSub[max_index] <= 2:
			exit("Error: too many processors {0:d} for lattice of size [{1:d},{1:d},{1:d},{2:d}]".format(numprocs, N, NT))
		if (restProc < 2):
			break


	return NSub[::-1]

# Dimensions
def get_processor_positions(N, NT, NP):
	max_lattice_size = N**3 * NT
	subdims = get_sub_dims(N, NT, NP)
	
	procs_per_dim = np.zeros(4)
	for i in xrange(3):
		procs_per_dim[i] = N / subdims[i]
	procs_per_dim[3] = NT / subdims[3]

	for i in xrange(3):
		assert procs_per_dim[i]*subdims[i] == N, "error in procs per dim"
	assert procs_per_dim[3]*subdims[3] == NT, "error in procs per dim"

	print "="*50
	print "N:                        ", N
	print "NT:                       ", NT
	print "Sub-dimensions:           ", subdims
	print "Processors per dimension: ", procs_per_dim
	print "="*50

	processor_positions = np.zeros((NP, 4), dtype=int)
	for rank in xrange(NP):
		processor_positions[rank][0] = rank % procs_per_dim[0]
		processor_positions[rank][1] = (rank / procs_per_dim[0]) % procs_per_dim[1]
		processor_positions[rank][2] = (rank / (procs_per_dim[0] * procs_per_dim[1])) % procs_per_dim[2]
		processor_positions[rank][3] = (rank / (procs_per_dim[0] * procs_per_dim[1] * procs_per_dim[2])) % procs_per_dim[3]

	for rank in xrange(NP):
		print "RANK: %3d POSITION: " % rank, processor_positions[rank]

	max_int_cpp = long(2147483647)
	SU3Doubles = 18 # number of su3 doubles
	SU3Size = 8 * SU3Doubles # number of bytes in su3
	linkDoubles = SU3Doubles * 4
	linkSize = linkDoubles * 8

	return processor_positions

def main():
	NSpatial = 48
	NTemporal = 96
	numprocs = 512

	get_processor_positions(NSpatial, NTemporal, numprocs)
	
if __name__ == '__main__':
	main()