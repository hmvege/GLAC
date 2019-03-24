import numpy as np
from operator import itemgetter

class color:
	# For easier printing
	PURPLE = '\033[95m'
	CYAN = '\033[96m'
	DARKCYAN = '\033[36m'
	BLUE = '\033[94m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
	END = '\033[0m'

class Index: 
	# General method for getting contiguous memory allocation
	def __init__(self,N):
		self.N = N
		self.N_length = len(N)
		
	def __call__(self,indices):
		offset = indices[0] + self.N[0]*(indices[1] + self.N[1]*(indices[2] + self.N[2]*indices[3]))
		return offset

class Processor:
	def __init__(self, rank):		
		self.rank = rank

	def set_neighbor_list(self, neighbor_list):
		if type(neighbor_list) != list:
			raise ValueError("Argument 'neighbor_list = %g' should be a list of processors." % neighbor_list)
		self.neighbor_list = neighbor_list

	def set_processor_position(self, proc_pos):
		if type(proc_pos) != list:
			raise ValueError("Argument 'proc_pos = %g' should be a list of processors." % proc_pos)
		self.proc_pos = proc_pos

	def __str__(self):
		string = "Rank %2d || x-1 = %2d x+1 = %2d || y-1 = %2d y+1 = %2d || z-1 = %2d z+1 = %2d || t-1 = %2d t+1 = %2d" % tuple([self.rank] + self.neighbor_list)
		return string


class Neighbors:
	def __init__(self, numprocs, n_proc_list):
		if type(n_proc_list) != list:
			raise ValueError("Argument 'numprocs = %g' should be a list of processors." % n_proc_list)
		self.n_proc_list = n_proc_list
		self.Np_x = n_proc_list[0]
		self.Np_y = n_proc_list[1]
		self.Np_z = n_proc_list[2]
		self.Np_t = n_proc_list[3]
		assert numprocs == self.Np_x*self.Np_y*self.Np_z*self.Np_t, "The number of processors %d do not match the lattice dimensions %s" % (numprocs," ".join([str(i) for i in n_proc_list]))
		self.processors = [Processor(rank) for rank in xrange(numprocs)]
		self.numprocs = numprocs

	def create_neighbor_lists(self):
		for proc in self.processors:
			n_list = [	self._getXMinusOne(proc.rank),
						self._getXPlusOne(proc.rank),
						self._getYMinusOne(proc.rank),
						self._getYPlusOne(proc.rank),
						self._getZMinusOne(proc.rank),
						self._getZPlusOne(proc.rank),
						self._getTMinusOne(proc.rank),
						self._getTPlusOne(proc.rank)]
			proc.set_neighbor_list(n_list)
			proc.set_processor_position(self._getProcessorIndexPosition(proc.rank))

	def print_basics(self):
		print "Number of processors: %3d" % self.numprocs
		print "Processors per dimensions: %s" % ' '.join([str(n) for n in self.n_proc_list])

	def print_neighbours(self,Np):
		n_list = self.processors[Np].neighbor_list
		print "Process rank = %d" % Np
		print "x-1 = %2d x+1 = %2d" % (n_list[0],n_list[1])
		print "y-1 = %2d y+1 = %2d" % (n_list[2],n_list[3])
		print "z-1 = %2d z+1 = %2d" % (n_list[4],n_list[5])
		print "t-1 = %2d t+1 = %2d" % (n_list[6],n_list[7])

	def print_processor_positions(self):
		print "\nPROCESSOR POSITIONS"
		for proc in self.processors:
			string = "|| %s%3d%s = (%2d, %2d, %2d, %2d)" % (color.BOLD,proc.rank,color.END,proc.proc_pos[0],proc.proc_pos[1],proc.proc_pos[2],proc.proc_pos[3])
			print string,
			self._print_break(proc.rank,string)

	def print_processor_lattice(self):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		for Np in xrange(numprocs):
			print "%2d" % Np,
			if (Np+1) % (Nx*Ny*Nz) == 0:
				print "\n"+6*"="
			elif (Np+1) % (Nx*Ny) == 0:
				print "\n"
			elif (Np+1) % Nx == 0:
				print ""

	def print_X_neighbors(self):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		print "\nX-DIRECTION"
		for Np in xrange(self.numprocs):
			# Method for +1
			neighbors_string = "|| %s%3d%s: x-1: %3d x+1: %3d" % (color.BOLD,Np,color.END,self.processors[Np].neighbor_list[0],self.processors[Np].neighbor_list[1])
			print neighbors_string,
			self._print_break(Np, neighbors_string)
		print "\n"

	def print_Y_neighbors(self):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		print "\nY-DIRECTION"
		for Np in xrange(self.numprocs):
			neighbors_string = "|| %s%3d%s: y-1: %3d y+1: %3d" % (color.BOLD,Np,color.END,self.processors[Np].neighbor_list[2],self.processors[Np].neighbor_list[3])
			print neighbors_string,
			self._print_break(Np, neighbors_string)
		print "\n"

	def print_Z_neighbors(self):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		print "\nZ-DIRECTION"
		for Np in xrange(self.numprocs):
			neighbors_string = "|| %s%3d%s: z-1: %3d z+1: %3d" % (color.BOLD,Np,color.END,self.processors[Np].neighbor_list[4],self.processors[Np].neighbor_list[5])
			print neighbors_string,
			self._print_break(Np, neighbors_string)
		print "\n"

	def print_T_neighbors(self):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		print "\nT-DIRECTION"
		for Np in xrange(self.numprocs):
			neighbors_string = "|| %s%3d%s: t-1: %3d t+1: %3d" % (color.BOLD,Np,color.END,self.processors[Np].neighbor_list[6],self.processors[Np].neighbor_list[7])
			print neighbors_string,
			self._print_break(Np, neighbors_string)
		print "\n"

	def _print_break(self, Np, neighbors_string):
		"""
		For printing breaks in the print-out of the neighboring processors.
		"""
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		if (Np+1) % (Nx*Ny*Nz) == 0:
			print "\n\n\n"
		elif (Np+1) % (Nx*Ny) == 0:
			print "\n"+"="*(len(neighbors_string)*Nx+3-Nx*len(color.BOLD+color.END))
		elif (Np+1) % Nx == 0:
			print ""

	def _getXPlusOne(self, Np):
		Nx = self.Np_x
		if ((Np + 1) % Nx) == 0:
			x_count = (Np/Nx * Nx)
			return x_count
		else:
			return Np+1

	def _getXMinusOne(self, Np):
		Nx = self.Np_x
		x_count = (Np-1) % Nx
		y_count = Np/Nx * Nx
		return x_count + y_count

	def _getYPlusOne(self, Np):
		Nx, Ny = self.Np_x, self.Np_y
		if ((Np/Nx + 1) % Ny) == 0:
			x_count = Np % Nx
			y_count = (Np / (Ny*Nx)) * (Nx*Ny) # Method for getting base layer by doing integer division
			return y_count + x_count
		else:
			return Np + Nx

	def _getYMinusOne(self, Np):
		Nx, Ny, Nz = self.Np_x, self.Np_y, self.Np_z
		x_count = Np % Nx
		y_count = (Np/Nx - 1) % Ny * Nx
		z_count = Np/(Nx*Ny) * (Nx*Ny)
		return x_count + y_count + z_count

	def _getZPlusOne(self, Np):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		if (((Np/(Nx*Ny) + 1) % Nz) == 0):
			x_count = Np % Nx
			y_count = (Np/Nx) % Ny * Nx
			z_count = 0
			t_count = Np/(Nx*Ny*Nz) % Nt * (Nx*Ny*Nz)
			return x_count + y_count + z_count + t_count
		else:
			return Np + Nx*Ny

	def _getZMinusOne(self, Np):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		x_count = Np % Nx
		y_count = Np / Nx % Ny * Nx
		z_count = (Np/(Nx*Ny) - 1) % Nz * (Nx*Ny)
		t_count = Np/(Nx*Ny*Nz) * (Nx*Ny*Nz)
		return x_count + y_count + z_count + t_count

	def _getTPlusOne(self, Np):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		if (((Np/(Nx*Ny*Nz) + 1) % Nt) == 0):
			x_count = Np % Nx
			y_count = (Np/Nx) % Ny * Nx
			z_count = (Np/(Nx*Ny)) % Nz * (Nx*Ny)
			return x_count + y_count + z_count
		else:
			return Np + Nx*Ny*Nz

	def _getTMinusOne(self, Np):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		x_count = Np % Nx
		y_count = ((Np / Nx) % Ny) * Nx
		z_count = (Np/(Nx*Ny) % Nz) * (Nx*Ny)
		t_count = (Np/(Nx*Ny*Nz) - 1) % Nt * (Nx*Ny*Nz)
		return x_count + y_count + z_count + t_count

	def _getProcessorIndexPosition(self,rank):
		Nx, Ny, Nz, Nt = self.Np_x, self.Np_y, self.Np_z, self.Np_t
		return [rank % Nx,
				(rank / Nx) % Ny,
				(rank / (Nx * Ny)) % Nz,
				(rank / (Nx * Ny * Nz)) % Nt]

def get_file_length(fname):
	# FOR COUNTING LINES IN FILE
	with open(fname,"r") as file:
		for i,l in enumerate(file):
			pass
	return i

def main():
	n_proc_list = [2,2,2,4]
	numprocs = np.prod(np.asarray(n_proc_list))
	neighbors = Neighbors(numprocs, n_proc_list)
	neighbors.print_basics()
	neighbors.create_neighbor_lists()

	for i in range(np.prod(n_proc_list)):
		neighbors.print_neighbours(i)
	# neighbors.print_processor_positions()

	# neighbors.print_X_neighbors()
	# neighbors.print_Y_neighbors()
	# neighbors.print_Z_neighbors()
	# neighbors.print_T_neighbors()

	# # Checking if indexes are correct
	# s = np.asarray(sorted(np.loadtxt("pros_pos_TEMP.txt",usecols=(1,3),dtype=str), key=lambda i: int(itemgetter(0)(i))))
	# s_positions = [[int(j) for j in i] for i in s[:,1]]
	# index_error = False
	# for proc, s_rank, s_pos in zip(neighbors.processors,map(int,s[:,0]),s_positions):
	# 	if proc.proc_pos != s_pos and proc.rank != s_rank:
	# 		print "ERROR: Theoretical: rank %2d = (%2d, %2d, %2d, %2d) C++ program: rank %2d = (%2d, %2d, %2d, %2d)" % (proc.rank, proc.proc_pos[0], proc.proc_pos[1], proc.proc_pos[2], proc.proc_pos[3], 
	# 																													s_rank, s_pos[0], s_pos[1], s_pos[2], s_pos[3])
	# 		index_error = True
	# if not index_error: print "INDEX TEST PASSED: Indexes matches theoretical values"

	# # Checking if we actually are accessing the correct neighbors
	# showErrorLine = False
	# raw_neighbors_file = "n_list_raw.txt"
	# file_length = get_file_length(raw_neighbors_file)
	# with open(raw_neighbors_file,"r") as raw_neighbor_lists:
	# 	for line,n in zip(raw_neighbor_lists,range(file_length)):
	# 		try:
	# 			rank, positive_index, positive_n, negative_index, negative_n = [int(i) for i in line.split()]
	# 			proc = neighbors.processors[rank]
	# 			if positive_n != proc.neighbor_list[positive_index] or negative_n != proc.neighbor_list[negative_index]:
	# 				print "%sLine %d%s C++: rank %2d positive(%2d): %2d negative(%2d): %d THEORY: %s" % (color.BOLD, n, color.END, rank,positive_index,positive_n,negative_index,negative_n,proc)#, "Theory: ", neighbors.processors[rank].neighbor_list 
	# 				if raw_input(""): continue
	# 			# print "Theory: ", neighbors.processors[rank].neighbor_list, "C++: ", rank, negative_n, positive_n
	# 		except ValueError:
	# 			if showErrorLine: print "Value error at line(%d): %s --> skipping line" % (n,line)
	# 			continue
	# 		except IndexError:
	# 			if showErrorLine: print "Index error at line(%d): %s --> skipping line" % (n,line)
	# 			continue


if __name__ == '__main__':
	main()