import numpy as np

class Index():
	def __init__(self,Nx,Ny,Nz):
		self.Nx = Nx
		self.Ny = Ny
		self.Nz = Nz

	def __call__(self,i,j,k):
		return self.Nz*(self.Ny*i+j)+k

def main():
	N = 3
	n = Index(N,N,N)
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



if __name__ == '__main__':
	main()