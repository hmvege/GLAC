import sympy as sp, numpy as np

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

def cMult(varName1, varName2, i,j,sign=False,conjugate=False):
	# z*w = (a + bi) * (c + di)
	# conjugate(z)*conjugate(w) = (a - bi) * (c - di) = ac - bd - i(ad + bc)
	# z*w = ac - bd + i(ad + bc)
	# - z*w =  - ac + bd - i(ad + bc)
	# conjugate: z*w = ac - bd - i(ad + bc)
	# conjugate: - z*w = - ac + bd + i(ad + bc)
	# (0  + 1i)  (2  + 3i)  (4  + 5i)
	# (6  + 7i)  (8  + 9i)  (10 + 11i)
	# (12 + 13i) (14 + 15i) (16 + 17i)
	a = 2*i
	b = 2*i+1
	c = 2*j
	d = 2*j+1
	if not conjugate:
		if not sign:
			msg_real = " + {0:<s}[{2:<d}]*{1:<s}[{4:<d}] - {0:<s}[{3:<d}]*{1:<s}[{5:<d}]".format(varName1,varName2,a,b,c,d)
			msg_imag = " + {0:<s}[{2:<d}]*{1:<s}[{5:<d}] + {0:<s}[{3:<d}]*{1:<s}[{4:<d}]".format(varName1,varName2,a,b,c,d)
		else:
			msg_real = " - {0:<s}[{2:<d}]*{1:<s}[{4:<d}] + {0:<s}[{3:<d}]*{1:<s}[{5:<d}]".format(varName1,varName2,a,b,c,d)
			msg_imag = " - {0:<s}[{2:<d}]*{1:<s}[{5:<d}] - {0:<s}[{3:<d}]*{1:<s}[{4:<d}]".format(varName1,varName2,a,b,c,d)
	else:
		if not sign:
			msg_real = " + {0:<s}[{2:<d}]*{1:<s}[{4:<d}] - {0:<s}[{3:<d}]*{1:<s}[{5:<d}]".format(varName1,varName2,a,b,c,d)
			msg_imag = " - {0:<s}[{2:<d}]*{1:<s}[{5:<d}] - {0:<s}[{3:<d}]*{1:<s}[{4:<d}]".format(varName1,varName2,a,b,c,d)
		else:
			msg_real = " - {0:<s}[{2:<d}]*{1:<s}[{4:<d}] + {0:<s}[{3:<d}]*{1:<s}[{5:<d}]".format(varName1,varName2,a,b,c,d)
			msg_imag = " + {0:<s}[{2:<d}]*{1:<s}[{5:<d}] + {0:<s}[{3:<d}]*{1:<s}[{4:<d}]".format(varName1,varName2,a,b,c,d)
	return msg_real,msg_imag

def printSympyCrossProductOld():
	for i in range(9):
		print "H%d," % i,
	print
	H0, H1, H2, H3, H4, H5, H6, H7, H8 = sp.symbols("H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8]")
	
	H = sp.Matrix([[H0, H1, H2], [H3, H4, H5], [H6, H7, H8]])
	H3Crossed = H[:,0].cross(H[:,1])
	for elem in H3Crossed:
		print elem

	exit()


def printSympyCrossProduct():
 	H0, H2, H4, H6, H8, H10, H12, H14, H16 = sp.symbols("H[0], H[2], H[4], H[6], H[8], H[10], H[12], H[14], H[16]")
 	H1, H3, H5, H7, H9, H11, H13, H15, H17 = sp.symbols("H[1], H[3], H[5], H[7], H[9], H[11], H[13], H[15], H[17]")
 	H1, H3, H5, H7, H9, H11, H13, H15, H17 = [-i*sp.I for i in [H1, H3, H5, H7, H9, H11, H13, H15, H17]]

 	H = sp.Matrix([	[H0 + H1, H2 + H3, H4 + H5],
 					[H6 + H7, H8 + H9, H10 + H11],
 					[H12 + H13, H14 + H15, H16 + H17]])


	H3Crossed = H[:,0].cross(H[:,1])
	H3Crossed = sp.expand(H3Crossed)
	
	print "//REAL PARTS OF CROSS PRODUCT: "
	for elem in H3Crossed:
		print sp.re(elem)
	print "//MAGINARY PARTS OF CROSS PRODUCT: "

	
	# print H[:,2]

	
def printCrossProduct():
	# H[2] = H[3].c()*H[7].c() - H[6].c()*H[4].c();  // Need complex multiplication here
 	# H[5] = H[6].c()*H[1].c() - H[0].c()*H[7].c();  // Need complex multiplication here
 	# H[8] = H[0].c()*H[4].c() - H[3].c()*H[1].c();  // Need complex multiplication here
 	c = True
 	# print cMult("A","B",0,0,False,c)
 	# for i in range(9):
 	# 	print "H[%d],"%(2*i+1),
 	# print 

 	H3H7 = cMult("H","H",4,6,conjugate=c)
 	H6H4 = cMult("H","H",7,3,True,conjugate=c)
 	H6H1 = cMult("H","H",7,0,conjugate=c)
 	H0H7 = cMult("H","H",1,6,True,conjugate=c)
 	H0H4 = cMult("H","H",1,3,conjugate=c)
 	H3H1 = cMult("H","H",4,0,True,conjugate=c)

 	H4 = H3H7[0] + H6H4[0]
 	H5 = H3H7[1] + H6H4[1]

 	H10 = H6H1[0] + H0H7[0]
 	H11 = H6H1[1] + H0H7[1]

 	H16 = H0H4[0] + H3H1[0]
 	H17 = H0H4[1] + H3H1[1]
 	print "H[4] =" + H4 + ";"
 	print "H[5] =" + H5 + ";"
 	print "H[10] =" + H10 + ";"
 	print "H[11] =" + H11 + ";"
 	print "H[16] =" + H16 + ";"
 	print "H[17] =" + H17 + ";"
 	


if __name__ == '__main__':
	# a*b = ac - bd + i(ad + bc)
	# printSympyCrossProductOld()
	printCrossProduct()
	# printRegularMatrix()
	# printContigiousComplexMatrixMultiplication(2)
	# printDeterminant()
	# printSympyCrossProduct()