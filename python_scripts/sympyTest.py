import numpy as np
import sympy as sp
import re

r11_r, r11_i, r12_r, r12_i, r21_r, r21_i, r22_r, r22_i = sp.symbols("r11re r11im r12re r12im r21re r21im r22re r22im")
s11_r, s11_i, s12_r, s12_i, s21_r, s21_i, s22_r, s22_i = sp.symbols("s11re s11im s12re s12im s21re s21im s22re s22im")
t11_r, t11_i, t12_r, t12_i, t21_r, t21_i, t22_r, t22_i = sp.symbols("t11re t11im t12re t12im t21re t21im t22re t22im")
r11_i, r12_i, r21_i, r22_i = [r*sp.I for r in [r11_i, r12_i, r21_i, r22_i]]
s11_i, s12_i, s21_i, s22_i = [s*sp.I for s in [s11_i, s12_i, s21_i, s22_i]]
t11_i, t12_i, t21_i, t22_i = [t*sp.I for t in [t11_i, t12_i, t21_i, t22_i]]


R = sp.Matrix([[r11_r + r11_i, r12_r + r12_i,0], [r21_r + r21_i, r22_r + r22_i,0], [0,0,1]])
S = sp.Matrix([[s11_r + s11_i, 0, s12_r + s12_i], [0,1,0], [s21_r + s21_i,0, s22_r + s22_i]])
T = sp.Matrix([[1,0,0],[0, t11_r + t11_i, t12_r + t12_i], [0, t21_r + t21_i, t22_r + t22_i]])


U = R*S*T
UExpanded = sp.expand(U)
# U_unsorted = [re.findall("([\+\-\ ]?)[^I]?([\w\*]{4,})",str(i)) for i in UExpanded]
U_unsorted = [re.findall("([\+\-])?[ ]?([\w\*]{4,})",str(i)) for i in UExpanded]

# print UExpanded[0]
# print U_unsorted[0]

# # Prints the regex values found
# for i, j in zip(range(len(UExpanded)),UExpanded):
# 	print "%d: " % i, j
# 	print "%d: " % i, re.findall("([\+\- ])[ ]?([\w\*]{5,})",str(j))

# # Prints the values of the U matrix, expanded
# for i,j in zip(range(len(sp.expand(U))),sp.expand(U)):
# 	print i
# 	sp.pprint(j)
# 	print "\n"

def expand_powers(line):
	return_vals = []
	for i in line:
		vals = [j for j in i.split("*") if j]
		if vals[-1] == "2":
			vals.pop(-1)
			vals.append(vals[-1])
			return_vals.append("*".join(vals))
		else:
			return_vals.append("*".join(vals))
	return return_vals

def convertToCpp(line):
	replaceDikt = { "11re": "[0]","11im": "[1]",
					"12re": "[2]","12im": "[3]",
					"21re": "[4]","21im": "[5]",
					"22re": "[6]","22im": "[7]"}
	return_values = []
	for elems in line:
		for key in replaceDikt.keys():
			if key in elems:
				elems = elems.replace(key,replaceDikt[key])
		return_values.append(elems)
	return " ".join(return_values)


reals = []
imags = []

for matrix_pos,i in zip(U_unsorted,range(len(U_unsorted))):
	line_imags = []
	line_reals = []
	for expr in matrix_pos:
		# print expr
		if expr[-1].split('*')[0] == "I":
			# print "*".join(expr[-1].split("*")[1:])
			line_imags.append(" ".join([expr[0],"*".join(expr[-1].split("*")[1:])]))
		else:
			line_reals.append(" ".join(expr))
	line_imags = expand_powers(line_imags)
	line_reals = expand_powers(line_reals)
	line_imags = convertToCpp(line_imags)
	line_reals = convertToCpp(line_reals)
	imags.append((2*i+1,line_imags))
	reals.append((2*i,line_reals))


for i,j in zip(reals, imags):
	print "X[%d] = %s;  \nX[%d] = %s;" % (i[0],"".join(i[1:][0]),j[0],"".join(j[1:][0]))
exit()

# sp.pprint(R)
# sp.pprint(S)
# sp.pprint(T)
# sp.pprint(U)

# r11_r, r11_i, r12_r, r12_i, r21_r, r21_i, r22_r, r22_i = sp.symbols("r_{11,re} r_{11,im} r_{12,re} r_{12,im} r_{21,re} r_{21,im} r_{22,re} r_{22,im}")
# s11_r, s11_i, s12_r, s12_i, s21_r, s21_i, s22_r, s22_i = sp.symbols("s_{11,re} s_{11,im} s_{12,re} s_{12,im} s_{21,re} s_{21,im} s_{22,re} s_{22,im}")
# t11_r, t11_i, t12_r, t12_i, t21_r, t21_i, t22_r, t22_i = sp.symbols("t_{11,re} t_{11,im} t_{12,re} t_{12,im} t_{21,re} t_{21,im} t_{22,re} t_{22,im}")

# r11, r12, r21, r22 = sp.symbols("r11 r12 r21 r22")
# s11, s12, s21, s22 = sp.symbols("s11 s12 s21 s22")
# t11, t12, t21, t22 = sp.symbols("t11 t12 t21 t22")

# R = sp.Matrix([[r11 ,r12,0],[r21,r22,0],[0,0,1]])
# S = sp.Matrix([[s11,0,s12],[0,1,0],[s21,0,s22]])
# T = sp.Matrix([[1,0,0],[0,t11,t12],[0,t21,t22]])

epsilon = 0.2
x = np.random.rand(4) - 0.5
x[0] = np.sign(x[0])*np.sqrt(1-epsilon*epsilon)
x[1:] = x[1:]*epsilon/np.linalg.norm(x[1:])
R = R.subs([(r11,np.complex(x[0],x[3])), (r12, np.complex(x[1],x[2])), (r21, np.complex(-x[1],x[2])), (r22, np.complex(x[0],-x[3]))])
U = U.subs([(r11,np.complex(x[0],x[3])), (r12, np.complex(x[1],x[2])), (r21, np.complex(-x[1],x[2])), (r22, np.complex(x[0],-x[3]))])

x = np.random.rand(4) - 0.5
x[0] = np.sign(x[0])*np.sqrt(1-epsilon*epsilon)
x[1:] = x[1:]*epsilon/np.linalg.norm(x[1:])
S = S.subs([(s11,np.complex(x[0],x[3])), (s12, np.complex(x[1],x[2])), (s21, np.complex(-x[1],x[2])), (s22, np.complex(x[0],-x[3]))])
U = U.subs([(s11,np.complex(x[0],x[3])), (s12, np.complex(x[1],x[2])), (s21, np.complex(-x[1],x[2])), (s22, np.complex(x[0],-x[3]))])

x = np.random.rand(4) - 0.5
x[0] = np.sign(x[0])*np.sqrt(1-epsilon*epsilon)
x[1:] = x[1:]*epsilon/np.linalg.norm(x[1:])
T = T.subs([(t11,np.complex(x[0],x[3])), (t12, np.complex(x[1],x[2])), (t21, np.complex(-x[1],x[2])), (t22, np.complex(x[0],-x[3]))])
U = U.subs([(t11,np.complex(x[0],x[3])), (t12, np.complex(x[1],x[2])), (t21, np.complex(-x[1],x[2])), (t22, np.complex(x[0],-x[3]))])


sp.pprint(sp.simplify(U.evalf()))
UH = U.H
sp.pprint(sp.simplify(UH*U.evalf()))

