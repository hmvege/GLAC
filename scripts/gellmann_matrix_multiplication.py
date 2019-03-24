import sympy as sp

# Setting Gell-Mann matrices
T = [sp.zeros(3,3) for i in range(8)]
# a=1
T[0][0,1] = 1
T[0][1,0] = 1
# a=2
T[1][0,1] = -sp.I
T[1][1,0] = sp.I
# a=3
T[2][0,0] = 1
T[2][1,1] = -1
# a=4
T[3][0,2] = 1
T[3][2,0] = 1
# a=5
T[4][0,2] = -sp.I
T[4][2,0] = sp.I
# a=6
T[5][1,2] = 1
T[5][2,1] = 1
# a=7
T[6][1,2] = -sp.I
T[6][2,1] = sp.I
# a=8
T[7][0,0] = 1
T[7][1,1] = 1
T[7][2,2] = -2
T[7] *= 1/sp.sqrt(3)

# Adding factor 1/2 to every generator
T = [Ta/2 for Ta in T]
T = [(-sp.I)*Ta for Ta in T]

for i,j in zip(T,T):
	sp.pprint((i*j).trace())
exit(1)

# Setting general matrix A
a11,a12,a13,a21,a22,a23,a31,a32,a33 = sp.symbols("a11 a12 a13 a21 a22 a23 a31 a32 a33",real=False)
A = sp.Matrix([[a11, a12, a13],[a21, a22, a23],[a31, a32, a33]])

R = [0 for i in range(8)]
for i,T_a in enumerate(T):
	R[i] = sp.re((T_a*A).trace())
	# sp.pprint((T_a*A).trace())
	# print R[i]
# exit(1)



res = sp.zeros(3,3)
for i,T_a in enumerate(T):
	res += T_a*R[i]

res = sp.simplify(res)

for i,elem in enumerate(res):
	print "\n\nELEMENT: %d " % (2*i)
	sp.pprint(sp.re(elem))
	print "\n\nELEMENT: %d " % (2*i+1)
	sp.pprint(sp.im(elem))

def printAll():
	# for i,T_a in enumerate(T):
	# 	print "T(%d):" % i 
	# 	sp.pprint(T_a)
	# sp.pprint(A)
	# sp.pprint(sp.simplify(res))
	sp.pprint(res)
	# sp.printing.print_ccode(res)
	# print "REAL:"
	# sp.pprint(sp.simplify(sp.re(-res*sp.I)))
	# print "IMAG:"
	# sp.pprint(sp.simplify(sp.im(-res*sp.I)))


index_map = """\nIndex map:
H =
0 1 2    0  1   2  3    4  5
3 4 5 =  6  7   8  9   10 11
6 7 8   12 13  14 15   16 17
	"""
print index_map

# printAll()