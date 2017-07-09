import numpy as np

U1 = np.matrix([[1+1j,1+2j,1+3j],[2+1j,2+2j,2+3j],[3+1j,3+2j,3+3j]])
U2 = np.matrix([[4+4j,4+5j,4+6j],[5+4j,4+5j,5+6j],[6+4j,6+5j,6+6j]])

print "Original matrices: \nU1 = \n", U1, "\nU2 = \n", U1

print "\nMultiplying: \n",U1*U2

print "\nAdding: \n", U1 + U2

print "\nSubtracting: \n", U1 - U2

print "\nTranspose of U1: \n", U1.T

print "\nConjugate of U1: \n", U1.conjugate()

print "\nComplex conjugate of U1: \n", U1.H