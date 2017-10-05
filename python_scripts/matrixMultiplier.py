def printRegularMatrix(N=3):
	msg = ""
	for i in range(N):
		for j in range(N):
			msg += "temp[%d] = " % (N*i+j)
			for k in range(N):
				msg += "mat[%d]*B.mat[%d]" % (N*i+k,N*k+j)
				if k != N-1: msg += " + "
			print msg + ";"
			msg = ""

#Trace:
#temp[0] = mat[0]*B.mat[0] + mat[1]*B.mat[3] + mat[2]*B.mat[6];
#temp[4] = mat[3]*B.mat[1] + mat[4]*B.mat[4] + mat[5]*B.mat[7];
#temp[8] = mat[6]*B.mat[2] + mat[7]*B.mat[5] + mat[8]*B.mat[8];

### Complex matrix multiplications
# 
# (0  + 1i)  (2  + 3i)  (4  + 5i)
# (6  + 7i)  (8  + 9i)  (10 + 11i)
# (12 + 13i) (14 + 15i) (16 + 17i)
# 0 1 2   11 12 13   00 01 02
# 3 4 5 = 21 22 23 = 10 11 12
# 6 7 8   31 32 33   20 21 22
#
# 0 1 = 11 12 = 00 01 = (0 + 1i) (2 + 3i)
# 2 3   21 22   10 11   (4 + 5i) (6 + 7i)
#
# zw = (a + bi)(c + id) = ac + iad + ibc - bd = ac - bd + i(ad + bc)

def complexContigiousMultiplication(N,i,j,k,verbose=False):
	a = 2*N*i + 2*k
	b = 2*N*i + 2*k + 1
	c = 2*N*k + 2*j
	d = 2*N*k + 2*j + 1
	msg_real = "mat[%d]*B[%d] - mat[%d]*B[%d]" % (a,c,b,d)
	msg_imag = "mat[%d]*B[%d] + mat[%d]*B[%d]" % (a,d,b,c)
	# msg_real = "A[%d]*B[%d] - A[%d]*B[%d]" % (a,c,b,d)
	# msg_imag = "A[%d]*B[%d] + A[%d]*B[%d]" % (a,d,b,c)
	# msg_real = "H[%d]*H[%d] - H[%d]*H[%d]" % (a,c,b,d)
	# msg_imag = "H[%d]*H[%d] + H[%d]*H[%d]" % (a,d,b,c)
	# msg_real = "U[%d]*U[%d] - U[%d]*U[%d]" % (a,c,b,d)
	# msg_imag = "U[%d]*U[%d] + U[%d]*U[%d]" % (a,d,b,c)
	if verbose: print "(Non-contigious) A[%d]*B[%d] = " % (i,j)
	return msg_real, msg_imag

def complexContigiousAddition(N,i,j,verbose=False):
	a = 2*i
	b = 2*i+1
	c = 2*j
	d = 2*j+1
	msg_real = "A[%d] + B[%d]" % (a,c)
	msg_imag = "A[%d] + B[%d]" % (b,d)
	if verbose: print "(Non-contigious) A[%d] + B[%d] = " % (i,j)
	return msg_real, msg_imag

def complexContigiousSubtraction(N,i,j,verbose=False):
	a = 2*i
	b = 2*i+1
	c = 2*j
	d = 2*j+1
	msg_real = "A[%d] - B[%d]" % (a,c)
	msg_imag = "A[%d] - B[%d]" % (b,d)
	if verbose: print "(Non-contigious) A[%d] - B[%d] = " % (i,j)
	return msg_real, msg_imag

def printContigiousComplexMatrixMultiplication(N=3):
	print "COMPLEX SU3 MATRIX: "
	msg = ""
	for i in range(N):
		for j in range(N):
			for k in range(2): # Complex variables
				msg += "temp[%d] = " % (2*i*N + 2*j + k)
				for l in range(N):
					msg += complexContigiousMultiplication(N,i,j,l)[k]
					if l != N-1: msg += " + ";
				msg += ";\n"
			print msg,
			msg = ""
	

def printDeterminant():
	N = 3
	z1 = complexContigiousMultiplication(N,0,0)
	z2 = complexContigiousMultiplication(N,1,3)
	z3 = complexContigiousMultiplication(N,2,6)
	w1 = ' + '.join([z1[0],z2[0],z3[0]])

	z4 = complexContigiousMultiplication(N,3,1)
	z5 = complexContigiousMultiplication(N,4,4)
	z6 = complexContigiousMultiplication(N,5,7)
	w2 = ' + '.join([z4[0],z5[0],z6[0]])

	z7 = complexContigiousMultiplication(N,6,2)
	z8 = complexContigiousMultiplication(N,7,5)
	z9 = complexContigiousMultiplication(N,8,8)
	w3 = ' + '.join([z7[0],z8[0],z9[0]])
	msg = w1 + ' + \n' + w2 + ' + \n' + w3
	print msg 


if __name__ == '__main__':
	# a*b = ac - bd + i(ad + bc)

	# printRegularMatrix()
	printContigiousComplexMatrixMultiplication(3)
	# printDeterminant()