import sympy as sp

x1,X12,X21,x2,X13,X31,x3,X23,X32 = sp.symbols("x1 X12 X21 x2 X13 X31 x3 X23 X32",real=False)
Y1 = sp.Matrix([[x1,X12,0],[X21,-x1,0],[0,0,0]])
Y2 = sp.Matrix([[x2,0,X13],[0,0,0],[X31,0,-x2]])
Y3 = sp.Matrix([[0,0,0],[0,x3,X23],[0,X32,-x3]])
Y1Inverse = (sp.eye(3)-0.25*Y1).inv()
Y2Inverse = (sp.eye(3)-0.25*Y2).inv()
Y3Inverse = (sp.eye(3)-0.5*Y3).inv()

N = 0
print "="*N
sp.pprint(sp.simplify(Y1Inverse))
print "="*N
sp.pprint(sp.simplify(Y2Inverse))
print "="*N
sp.pprint(sp.simplify(Y3Inverse))

index_map = """
Index map:
H =
0 1 2    0  1   2  3    4  5
3 4 5 =  6  7   8  9   10 11
6 7 8   12 13  14 15   16 17
"""
print index_map
