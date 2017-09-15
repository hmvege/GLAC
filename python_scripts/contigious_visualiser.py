import numpy as np, sys

class Index: # General method for getting contiguous memory allocation
	def __init__(self,N):
		self.N = N
		self.N_length = len(N)
		
	def __call__(self,indices):
		if len(indices) != self.N_length:
			raise ValueError("Number of indices %d do not correspond with number of dimensions %d." % (len(indices),self.N_length))
		offset = int(0);
		for i in range(self.N_length):
			prod = 1
			for j in range(i+1,self.N_length):
				prod *= self.N[j]
			offset += indices[i]*prod
		return offset

class IndexOld: # Old method
	def __init__(self,Nx,Ny,Nz):
		self.Nx = Nx
		self.Ny = Ny
		self.Nz = Nz

	def __call__(self,i,j,k):
		return self.Nz*(self.Ny*i+j)+k

class Lattice:
	def __init__(self,N):
		if type(N) != list: raise ValueError("Dimension to be passed as elements in a list.")
		self.N = N

	def printLattice(self,reversed=False):
		# Redirecting output
		orig_stdout = sys.stdout
		f = open('lattice.txt', 'w')
		sys.stdout = f

		n = Index(self.N)
		memory_position_string = []
		if not reversed:
			msg = "XYZT order"
		else:
			msg = "TZYX order"
		scr_len = len("[%d,%d,%d,%d = %3d]" % (0,0,0,0,n([0,0,0,0])))*4+4
		str_len = ((len("[%d,%d,%d,%d = %3d]" % (0,0,0,0,n([0,0,0,0])))*4+3) - len(msg))/2
		print "="*str_len + msg + "="*(str_len+1)

		if not reversed:
			for x in xrange(self.N[0]):
				print "X = %d," % x,
				for y in xrange(self.N[1]):
					print "Y = %d" % y
					for z in xrange(self.N[2]):
						for t in xrange(self.N[3]):
							memory_position_string.append(n([x,y,z,t]))
							print_string = "[%d,%d,%d,%d = %3d]" % (x,y,z,t,memory_position_string[-1])
							print print_string,
						if ((x+1)%(self.N[0]) == 0 and (y+1)%(self.N[1]) == 0 and (z+1)%(self.N[2]) == 0 and (t+1)%(self.N[3]) == 0):
							print "\n","="*(scr_len-1)
						print "\n",
					print "\n",
				print "\n",
		else:
			for t in xrange(self.N[3]):
				print "T = %d," % t,
				for z in xrange(self.N[2]):
					print "Z = %d" % z
					for y in xrange(self.N[1]):
						for x in xrange(self.N[0]):
							memory_position_string.append(n([x,y,z,t]))
							print_string = "[%d,%d,%d,%d = %3d]" % (x,y,z,t,memory_position_string[-1])
							print print_string,
						if ((x+1)%(self.N[0]) == 0 and (y+1)%(self.N[1]) == 0 and (z+1)%(self.N[2]) == 0 and (t+1)%(self.N[3]) == 0):
							print "\n","="*(scr_len-1)
						print "\n",
					print "\n",
				print "\n",
		
		# Closing output file
		sys.stdout = orig_stdout
		f.close()

		# Storing memory string
		self.memory_position_string = memory_position_string



class SubLattice:
	def __init__(self,N,numprocs):
		if type(N) != list: raise ValueError("Dimension to be passed as elements in a list.")
		self.N = [n for n in N]
		self.NSub = [n for n in N]
		self.numprocs = numprocs
		restProc = numprocs
		self.nSub = Index(self.NSub)
		self.n = Index(self.N)

		# Sets up the sub-lattice dimensions
		while restProc >= 2:
			for i in range(0,4):
				self.NSub[i] /= 2;
				restProc /= 2;
				if (restProc < 2): break

		# Sets up processors per dimension
		self.procsPerDim = []
		for i in range(4):
			self.procsPerDim.append(self.N[i] / self.NSub[i])

		# Sets up the volumes
		self.V, self.VSub, self.VProc, = [],[],[]
		self.V.append(self.N[0])
		self.VSub.append(self.NSub[0])
		self.VProc.append(self.procsPerDim[0])
		for i in range(3):
			self.V.append(self.N[i]*self.V[-1])
			self.VSub.append(self.NSub[i]*self.VSub[-1])
			self.VProc.append(self.procsPerDim[i]*self.VProc[-1])

		# Map for 2^4 prosessors
		# t=0, z = 0
		#  0  1 # y = 0
		#  2  3 # y = 1
		# z = 1
		#  4  5 # y = 0
		#  6  7 # y = 1

		# t = 1, z = 0
		#  8  9 # y = 0
		# 10 11 # y = 1
		# z = 1
		# 12 13 # y = 0
		# 14 15 # y = 1

		# Sets up processor positions in the lattice
		self.processor_coordinates = []
		for rank in range(numprocs):
			processor_coordinate = [rank % self.procsPerDim[0],
									(rank / self.VProc[0]) % self.procsPerDim[1],
									(rank / self.VProc[1]) % self.procsPerDim[2],
									(rank / self.VProc[2]) % self.procsPerDim[3]]
			self.processor_coordinates.append(processor_coordinate)

		print "Lattice dimensions: ", self.N
		print "Volumes of lattice: ", self.V
		print "Sub-lattice dimensions: ", self.NSub
		print "Volumes of sub-lattice: ", self.VSub
		print "Number of processsors: ", self.numprocs
		print "Processors per dimension: ", self.procsPerDim
		print "Volumes of processsors: ", self.VProc
		# print "Processor coordinates: \n",self.processor_coordinates

	def printTotalLattice(self):
		# Redirecting output
		orig_stdout = sys.stdout
		f = open('all_sub_lattices.txt', 'w')
		sys.stdout = f

		memory_position_string = [] # This can be compared later with the scalar version. If they are equal, the methods are working and the bug is not due to wrong memory placement.
		x_offset, y_offset, z_offset, t_offset = 0,0,0,0

		for prosessor in range(self.numprocs):
			str_processor =  "Prosessor: %2d" % prosessor
			print str_processor
			msg = "TZYX order"
			scr_len = len("[%d,%d,%d,%d = %3d]" % (0,0,0,0,self.nSub([1,1,1,1])))*self.NSub[0]+self.NSub[0]
			str_len = ((len("[%d,%d,%d,%d = %3d]" % (0,0,0,0,self.nSub([1,1,1,1])))*self.NSub[0]+1) - len(msg))/2
			str_header = "="*str_len + msg + "="*(str_len+1)
			print str_header

# nt = m_V[2] * (m_neighbourLists->getProcessorDimensionPosition(3) * m_VSub[3] + t) * linkSize;
#         for (int z = 0; z < m_NSpatial; z++) {
#             nz = m_V[1] * (m_neighbourLists->getProcessorDimensionPosition(2) * m_VSub[2] + z) * linkSize + nt;
#             for (int y = 0; y < m_NSpatial; y++) {
#                 ny = m_V[0] * (m_neighbourLists->getProcessorDimensionPosition(1) * m_VSub[1] + y) * linkSize + nz;
#                 for (int x = 0; x < m_NSpatial; x++) {
#                     nx = (m_neighbourLists->getProcessorDimensionPosition(0) * m_VSub[0] + x) * linkSize + ny;
#                     startPoints = nx;
#                     MPI_File_write_at(file,startPoints, &m_lattice[m_indexHandler->getIndex(x,y,z,t)], 72*sizeof(double), MPI_DOUBLE, MPI_STATUS_IGNORE);
#                 }
#             }
#         }

			for t in xrange(self.NSub[3]):
				t_offset = self.V[2] * (self.processor_coordinates[prosessor][3] * self.NSub[3] + t)
				str_t_cord = "T = %d," % t
				print str_t_cord,
				for z in xrange(self.NSub[2]):
					z_offset = self.V[1] * (self.processor_coordinates[prosessor][2] * self.NSub[2] + z) + t_offset
					str_z_cord = "Z = %d" % z
					print str_z_cord
					for y in xrange(self.NSub[1]):
						y_offset = self.V[0] * (self.processor_coordinates[prosessor][1] * self.NSub[1] + y) + z_offset
						for x in xrange(self.NSub[0]):
							x_offset = self.processor_coordinates[prosessor][0] * self.NSub[0] + x + y_offset
							print_string = "[%d,%d,%d,%d = %3d]" % (x,y,z,t,x_offset) # Memory positions
							print print_string,
						if ((x+1)%(self.NSub[0]) == 0 and (y+1)%(self.NSub[1]) == 0 and (z+1)%(self.NSub[2]) == 0 and (t+1)%(self.NSub[3]) == 0):
							print "\n","="*(scr_len-1)
						print "\n",
					print "\n",
				print "\n",

		self.memory_position_string = memory_position_string
		sys.stdout = orig_stdout
		f.close()


	def printSubLattice(self, reversed = False):
		# Redirecting output
		orig_stdout = sys.stdout
		f = open('sublattice.txt', 'w')
		sys.stdout = f

		for prosessor in range(self.numprocs):
			print "Prosessor: ", prosessor
			if not reversed:
				msg = "XYZT order"
			else:
				msg = "TZYX order"
			scr_len = len("[%d,%d,%d,%d = %3d]" % (0,0,0,0,self.nSub([0,0,0,0])))*self.NSub[0]+self.NSub[0]
			str_len = ((len("[%d,%d,%d,%d = %3d]" % (0,0,0,0,self.nSub([0,0,0,0])))*self.NSub[0]+1) - len(msg))/2
			print "="*str_len + msg + "="*(str_len+1)

			if not reversed:
				for x in xrange(self.NSub[0]):
					print "X = %d," % x,
					for y in xrange(self.NSub[1]):
						print "Y = %d" % y
						for z in xrange(self.NSub[2]):
							for t in xrange(self.NSub[3]):
								print_string = "[%d,%d,%d,%d = %3d]" % (x,y,z,t,self.nSub([x,y,z,t]))
								print print_string,
							if ((x+1)%(self.NSub[0]) == 0 and (y+1)%(self.NSub[1]) == 0 and (z+1)%(self.NSub[2]) == 0 and (t+1)%(self.NSub[3]) == 0):
								print "\n","="*(scr_len-1)
							print "\n",
						print "\n",
					print "\n",
			else:
				for t in xrange(self.NSub[3]):
					print "T = %d," % t,
					for z in xrange(self.NSub[2]):
						print "Z = %d" % z
						for y in xrange(self.NSub[1]):
							for x in xrange(self.NSub[0]):
								print_string = "[%d,%d,%d,%d = %3d]" % (x,y,z,t,self.nSub([x,y,z,t]))
								print print_string,
							if ((x+1)%(self.NSub[0]) == 0 and (y+1)%(self.NSub[1]) == 0 and (z+1)%(self.NSub[2]) == 0 and (t+1)%(self.NSub[3]) == 0):
								print "\n","="*(scr_len-1)
							print "\n",
						print "\n",
					print "\n",

		sys.stdout = orig_stdout
		f.close()


def old_index_test(N=3):
	print "Old index test"
	n = IndexOld(N,N,N)
	M = np.zeros((N,N,N,3),dtype=tuple)
	for i in xrange(N):
		for j in xrange(N):
			for k in xrange(N):
				M[i,j,k] = (i,j,k)
				# print n(i,j,k)
	
	for i in xrange(N): # X
		for j in xrange(N): # Y
			print "[",
			for k in xrange(N): # Z
				print "%d %d %d: n=%2d" % (M[i,j,k][0],M[i,j,k][1],M[i,j,k][2],n(i,j,k)),
				if k < N-1: print ",",
			print "]\n",
		print ""

def new_index_test(N):
	if type(N) != list: raise TypeError("Please provide a list of dimension sizes.")
	print "New index test"
	n = Index(N)
	M = np.zeros((N[0],N[1],N[2],3),dtype=tuple)
	for i in xrange(N[0]):
		for j in xrange(N[1]):
			for k in xrange(N[2]):
				M[i,j,k] = (i,j,k)
				# print n(i,j,k)
	
	for i in xrange(N[0]): # X
		for j in xrange(N[1]): # Y
			print "[",
			for k in xrange(N[2]): # Z
				print "%d %d %d: n=%2d" % (M[i,j,k][0],M[i,j,k][1],M[i,j,k][2],n([i,j,k])),
				if k < N[2]-1: print ",",
			print "]\n",
		print ""


def main():
	# old_index_test()
	# new_index_test([3,3,3])
	dimensions = [4,4,4,4]
	numprocs = 16
	lattice = Lattice(dimensions)
	lattice.printLattice(True)
	sublattice = SubLattice(dimensions,numprocs)
	# sublattice.printSubLattice()
	sublattice.printTotalLattice()
	if sublattice.memory_position_string == lattice.memory_position_string:
		print "Success"
	else:
		print "Not equal memory positions"

if __name__ == '__main__':
	main()