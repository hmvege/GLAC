import sympy as sp



x1,X12,X21,x2,X13,X31,x3,X23,X32 = sp.symbols("x1 X12 X21 x2 X13 X31 x3 X23 X32",real=False)
Y1 = sp.Matrix([[x1,X12,0],[X21,-x1,0],[0,0,0]])
Y2 = sp.Matrix([[x2,0,X13],[0,0,0],[X31,0,-x2]])
Y3 = sp.Matrix([[0,0,0],[0,x3,X23],[0,X32,-x3]])
I = sp.eye(3)

# Finding the inverses
Y1Inverse = (1.0*I - 0.25*Y1).inv()
Y2Inverse = (1.0*I - 0.25*Y2).inv()
Y3Inverse = (1.0*I - 0.5*Y3).inv()

# Setting up the U's
U1 = (1.0*I + 0.25*Y1)*Y1Inverse
U2 = (1.0*I + 0.25*Y2)*Y2Inverse
U3 = (1.0*I + 0.5*Y3)*Y3Inverse

# Simplifying the U's
U1 = sp.simplify(U1)
U2 = sp.simplify(U2)
U3 = sp.simplify(U3)

# def printSeperator(length=100): print length*"="

# for i in range(3):
# 	for j in range(3):
# 		printSeperator()
# 		print "\nIndex: %d" % (6*i + 2*j)
# 		sp.pprint(sp.simplify(sp.re(U1[i,j])))
# 		printSeperator()
# 		print "\nIndex: %d" % (6*i + 2*j + 1)
# 		sp.pprint(sp.simplify(sp.im(U1[i,j])))

# Printing results
print "\nU1: "
sp.pprint(U1)
print "\nU2: "
sp.pprint(U2)
print "\nU3: "
sp.pprint(U3)



# sp.pprint(sp.simplify(Y1Inverse))
# print "\n"
# sp.pprint(sp.simplify(Y2Inverse))
# print "\n"
# sp.pprint(sp.simplify(Y3Inverse))

index_map = """
Index map:
H =
0 1 2    0  1   2  3    4  5
3 4 5 =  6  7   8  9   10 11
6 7 8   12 13  14 15   16 17
"""
print index_map
