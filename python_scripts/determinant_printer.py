import sympy as sp
import numpy as np
from sympy.matrices.expressions import Determinant
import re

# H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15, H16, H17 = sp.symbols("""H0 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17""")
# H1, H3, H5, H7, H9, H11, H13, H15, H17 = [i*sp.I for i in [H1, H3, H5, H7, H9, H11, H13, H15, H17]]

H0, H1, H2, H3, H4, H5, H6, H7 = sp.symbols("H0, H1, H2, H3, H4, H5, H6, H7")
H1, H3, H5, H7 = [i*sp.I for i in [H1, H3, H5, H7]]

def print_matrix_line(N=3,print_complex=False):
	msg = ""
	for i in range(N):
		for j in range(N):
			if print_complex:
				msg += ("H%d, " % (2*N*i + 2*j + 1))
			else:
				for k in range(2):
					msg += ("H%d, " % (2*N*i + 2*j + k))
	print msg

def print_matrix(N=3):
	msg = "["
	for i in range(N):
		msg += "["
		for j in range(N):
			for k in range(2):
				msg += "H%d" % (2*N*i + 2*j + k)
				if k != 1: msg += "+"
			if j != (N-1): msg += ","
		msg += "]"
		if i != (N-1): msg += ","
	msg += "]"
	print msg

def convertToCpp(string):
	replaceDikt = {	"0": "[0]","1": "[1]","2": "[2]","3": "[3]","4": "[4]","5": "[5]","6": "[6]","7": "[7]","8": "[8]",
					"9": "[9]","10": "[10]","11": "[11]","12": "[12]","13": "[13]","14": "[14]","15": "[15]","16": "[16]","17": "[17]"}
	return_values = []
	string = string.split(" ")
	string = [s.split("*") for s in string]
	for elems in string:
		for atom in elems:
			for key in replaceDikt.keys():
				if key in atom:
					# print key, atom
					atom = atom.replace(key,replaceDikt[key])
			return_values.append(atom)

	# return " ".join([" * ".join(s) for s in string])
	return " ".join(return_values)

# Small function for writing out dictionaries in text that can be copied into python as code.
# s = "{"
# for i in range(18):
# 	s += "\"%d\": \"[%d]\"" % (i,i)
# 	if i != 17: s += ","
# s += "}"
# print s

# print_matrix_line(2)
# print_matrix(2)

# H = sp.Matrix([[H0+H1,H2+H3,H4+H5],[H6+H7,H8+H9,H10+H11],[H12+H13,H14+H15,H16+H17]])
H = sp.Matrix([[H0+H1,H2+H3],[H4+H5,H6+H7]])
HDet = H.det()



# print HDet
H_unsorted = re.findall("([\+\-])?[ ]?([\w\*]{4,})",str(HDet))
H_real = []
H_imag = []
for H_element in H_unsorted:
	if H_element[-1].split("*")[0] == "I":
		H_imag.append(" ".join([H_element[0],"*".join(H_element[-1].split("*")[1:])]))
	else:
		H_real.append(" ".join(H_element))

print "det(%s,\n%s);" % (" ".join(H_real)," ".join(H_imag))
H_str_real = convertToCpp(" ".join(H_real))
H_str_imag = convertToCpp(" ".join(H_imag))

print "det(%s,\n%s);" % (H_str_real,H_str_imag)



