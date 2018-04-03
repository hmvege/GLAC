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
	print "Trace: ", np.trace(U1)

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

	print "// Hermitian U1 matrix"
	print write_to_cpp((U1+U1.H) - np.trace((U1+U1.H))*0.5,"U1Hermitian",3)

	print "// Anti-hermitian U1 matrix"
	print write_to_cpp((U1-U1.H) - np.trace((U1-U1.H))*0.5,"U1AntiHermitian",3)


def write_to_cpp(mat,mat_name,N):
	string = ""
	for i in xrange(N):
		for j in xrange(N):
			string += "%s.setComplex(complex(%g,%g),%d);" % (mat_name,np.real(mat[i,j]),np.imag(mat[i,j]),int(2*N*i + 2*j))
			# if i != 1 and j != 1: string += "\n"
			if (j!=N-1): string += "\n"
		if (i!=N-1): string += "\n"
	# exit(1)
	return string

def write_complex_to_cpp(w,name):
	string = ""
	if (type(w)==complex):
		string = "z%s = complex(%g,%g);" % (name,w.real,w.imag)
	elif (type(w)==np.float64 or type(w)==float):
		string = "z%s = %.15g;" % (name,w)
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

def complex_class_unit_tests():
	z1 = 1.0 + 2.0j
	z2 = 3.0 + 4.0j

	print "// Adding"
	print write_complex_to_cpp(z1+z2,"Add")
	print "// Subtracting"
	print write_complex_to_cpp(z1-z2,"Sub")
	print "// Multiplying"
	print write_complex_to_cpp(z1*z2,"Mul")
	print "// Division"
	print write_complex_to_cpp(z1/z2,"Div")
	print "// Conjugate 1, conjugate()/c()"
	print write_complex_to_cpp(z1.conjugate(),"Conj")
	# complex.
	# print "// Conjugate 2, c()"
	# print write_complex_to_cpp(z1.conjugate())
	print "// Norm"
	print write_complex_to_cpp(np.linalg.norm(z1),"Norm")
	print "// Norm squared"
	print write_complex_to_cpp(np.linalg.norm(z1)**2,"NormSquared")
	print "// Set to minus operator"
	print write_complex_to_cpp(-z1,"SetToMinus")

if __name__ == '__main__':
	write_old_matrices()
	write_matrices_contigious_complex_su3()
	complex_class_unit_tests()
