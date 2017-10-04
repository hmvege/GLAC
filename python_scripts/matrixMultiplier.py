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
# (0  + 1i)  (2  + 3i)  (4  + 5i)
# (6  + 7i)  (8  + 9i)  (10 + 11i)
# (12 + 13i) (14 + 15i) (16 + 17i)
#
# 0 1 2   11 12 13   00 01 02
# 3 4 5 = 21 22 23 = 10 11 12
# 6 7 8   31 32 33   20 21 22
# 

def printContigiousComplexMatrixMultiplication(N=3):
	msg = ""
	for i in range(N):
		for j in range(N):
			# print "M[%d,%d] = %2d |" % (i,j,2*i*N+2*j),
			msg += "%d: " % (N*i + j)
			for k in range(2): # Complex variables
				msg += "%3d" % (2*i*N + 2*j + k)
			if j != (N-1): msg += " | "
		print msg,""
		msg = ""

if __name__ == '__main__':
	# printRegularMatrix()
	printContigiousComplexMatrixMultiplication()
