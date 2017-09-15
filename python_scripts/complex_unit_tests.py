import numpy as np

def write_to_cpp(mat,mat_name):
	string = ""
	for i in xrange(3):
		for j in xrange(3):
			string += "%s.mat[%d] = complex(%g,%g);\n" % (mat_name,int(i*3+j),np.real(mat[i,j]),np.imag(mat[i,j]))
	return string

U1 = np.matrix([[1+1j,1+2j,1+3j],[2+1j,2+2j,2+3j],[3+1j,3+2j,3+3j]])
U2 = np.matrix([[4+4j,4+5j,4+6j],[5+4j,5+5j,5+6j],[6+4j,6+5j,6+6j]])

print "Original matrices: \nU1 = \n", U1, "\nU2 = \n", U1

print "\nAdding:"
print write_to_cpp(U1+U2,"UAdd")

print "\nSubtracting:"
print write_to_cpp(U1-U2,"USub")

print "\nMultiplying:"
print write_to_cpp(U1*U2,"UMul")

print "\nConjugate of U1:"
print write_to_cpp(U1.conjugate(),"UC")

print "\nTranspose of U1:"
print write_to_cpp(U1.T,"UT")

print "\nComplex conjugate of U1:"
print write_to_cpp(U1.H,"UCT")
