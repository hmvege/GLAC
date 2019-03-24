import numpy as np

U1 = np.asarray([[   1 +    1j ,     1 +    2j ,     1 +    3j   ],
[   2 +    1j ,     2 +    2j ,     2 +    3j   ],
[   3 +    1j ,     3 +    2j ,     3 +    3j   ]],dtype=complex)

U2 = np.asarray([[   4 +    4j ,     4 +    5j ,     4 +    6j  ],
[   5 +    4j ,     5 +    5j ,     5 +    6j  ],
[   6 +    4j ,     6 +    5j ,     6 +    6j  ]],dtype=complex)

U3 = np.asarray([[   4 +    4j ,     5 +    5j ,     2 +    6j  ],
[   6 +    4j ,     5 +    5j ,     1 +    6j  ],
[   6 +    4j ,     9 +    5j ,     2 +    6j  ]],dtype=complex)

UAdd = U1 + U2

r = np.asarray([[   1 +    2j ,     3 +    4j],
				[   5 +    6j ,     7 +    8j]],dtype=complex)
s = np.asarray([[   2 +    5j ,     4 +    8j],
 				[   2 +    6j ,  	10+    3j]],dtype=complex)
t = np.asarray([[   7 +    2j ,     6 +    1j],
				[   6 +    4j ,     5 +    2j]],dtype=complex)

R = np.asarray([[r[0,0],r[0,1],0],
				[r[1,0],r[1,1],0],
				[0,0,1]],dtype=complex)

S = np.asarray([[s[0,0],0,s[0,1]],
				[0,1,0],
				[s[1,0],0,s[1,1]]],dtype=complex)

T = np.asarray([[1,0,0],
				[0,t[0,0],t[0,1]],
				[0,t[1,0],t[1,1]]],dtype=complex)

def printMatrix(A,symbol = "H", N=3):
	for i in range(N):
		for j in range(N):
			print symbol + ".setComplex(complex(%d,%d),%d);" % (np.real(A[i,j]),np.imag(A[i,j]), 2*N*i + 2*j) 
			# print symbol + "[%d] = " % (2*N*i + 2*j) + "%g" % np.real(A[i,j]) + ";"
			# print symbol + "[%d] = " % (2*N*i + 2*j + 1) + "%g" % np.imag(A[i,j]) + ";"

def printDotProduct():
	print U1,'\n'

	print U2,'\n'

	print np.dot(U1,U2)

if __name__ == '__main__':
	# printDotProduct()
	U_traced = np.dot(U1,U3)
	RST = np.dot(np.dot(R,S),T)
	printMatrix(r,"s_r",2)
	printMatrix(s,"s_s",2)
	printMatrix(t,"s_t",2)
	printMatrix(R,"U_R",3)
	printMatrix(S,"U_S",3)
	printMatrix(T,"U_T",3)
	printMatrix(RST,"URST",3)


	# print U_traced
	# print np.matrix.trace(U_traced)

	# printMatrix(U1,"U1")
	# printMatrix(U2,"U2")
	# printMatrix(U3,"UTrace")
	# printMatrix(U1*U2*UAdd,"URST")