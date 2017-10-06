import numpy as np

# OLD PROGRAM
def write_to_cpp_old(mat,mat_name):
	string = ""
	for i in xrange(3):
		for j in xrange(3):
			string += "%s.mat[%d] = complex(%g,%g);\n" % (mat_name,int(i*3+j),np.real(mat[i,j]),np.imag(mat[i,j]))
	return string

def write_old_matrices():
	U1 = np.matrix([[1+1j,1+2j,1+3j],[2+1j,2+2j,2+3j],[3+1j,3+2j,3+3j]])
	U2 = np.matrix([[4+4j,4+5j,4+6j],[5+4j,5+5j,5+6j],[6+4j,6+5j,6+6j]])

	print "Original matrices: \nU1 = \n", U1, "\nU2 = \n", U2

	print "// Adding"
	print write_to_cpp_old(U1+U2,"UAdd")

	print "// Subtracting"
	print write_to_cpp_old(U1-U2,"USub")

	print "// Multiplying"
	print write_to_cpp_old(U1*U2,"UMul")

	print "// Conjugate of U1"
	print write_to_cpp_old(U1.conjugate(),"UC")

	print "// Transpose of U1"
	print write_to_cpp_old(U1.T,"UT")

	print "// Complex conjugate of U1"
	print write_to_cpp_old(U1.H,"UCT")

def write_to_cpp(mat,mat_name,N):
	string = ""
	for i in xrange(N):
		for j in xrange(N):
			string += "%s.setComplex(complex(%g,%g),%d);" % (mat_name,np.real(mat[i,j]),np.imag(mat[i,j]),int(2*N*i + 2*j))
			# if i != 1 and j != 1: string += "\n"
			if (j!=1): string += "\n"
		if (i!=1): string += "\n"
	# exit(1)
	return string

def write_matrices_contigious_complex_su3():
	s1 = np.matrix([[1+1j,1+2j],[2+1j,2+2j]])
	s2 = np.matrix([[4+4j,4+5j],[5+4j,5+5j]])
	N = 2
	print "Original matrices:"
	print "// s1"
	print write_to_cpp(s1,"s1",N)

	print "// s2"
	print write_to_cpp(s2,"s2",N)

	print "// Adding"
	print write_to_cpp(s1+s2,"sAdd",N)

	print "// Subtracting"
	print write_to_cpp(s1-s2,"sSub",N)

	print "// Multiplying"
	print write_to_cpp(s1*s2,"sMul",N)

	print "// Conjugate of s1"
	print write_to_cpp(s1.conjugate(),"sC",N)

	print "// Transpose of s1"
	print write_to_cpp(s1.T,"sT",N)

	print "// Complex conjugate of s1"
	print write_to_cpp(s1.H,"sCT",N)

if __name__ == '__main__':
	# write_old_matrices()
	write_matrices_contigious_complex_su3()