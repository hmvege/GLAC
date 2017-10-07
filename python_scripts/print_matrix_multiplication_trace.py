import numpy as np
import sympy as sp
import re

def print_matrix_line(N=3,symbol="H",print_complex=False):
	msg = ""
	for i in range(N):
		for j in range(N):
			if print_complex:
				msg += ("%s%d, " % (symbol,2*N*i + 2*j + 1))
			else:
				for k in range(2):
					msg += ("%s%d, " % (symbol,2*N*i + 2*j + k))
	print msg

def print_matrix(N=3,symbol="H",):
	msg = "["
	for i in range(N):
		msg += "["
		for j in range(N):
			for k in range(2):
				msg += "%s%d" % (symbol,2*N*i + 2*j + k)
				if k != 1: msg += "+"
			if j != (N-1): msg += ","
		msg += "]"
		if i != (N-1): msg += ","
	msg += "]"
	print msg

def convertToCpp(line,N=3,symbol=["H"]):
	# Creating the dictionary - danger, eval!
	raw_dict = "{"
	for sym,j in zip(symbol,range(len(symbol))):
		for i in range(2*N*N):
			raw_dict += "\"{0:<s}{1:<d}\": \"{0:<s}[{1:<d}]\"".format(sym,i)
			if (i != (2*N*N - 1)) or (j != (len(symbol)-1)):
				raw_dict += ","
	raw_dict += "}"

	replaceDict = eval(raw_dict)
	return_values = []
	for elems in line:
		for key in replaceDict.keys():
			if key in elems:
				elems = elems.replace(key,replaceDict[key])
		return_values.append(elems)
	return "".join(return_values)

def get_trace_matrix():
	A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17 = sp.symbols("A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17")
	B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17 = sp.symbols("B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17")
	A1, A3, A5, A7, A9, A11, A13, A15, A17 = [i*sp.I for i in [A1, A3, A5, A7, A9, A11, A13, A15, A17]]
	B1, B3, B5, B7, B9, B11, B13, B15, B17 = [i*sp.I for i in [B1, B3, B5, B7, B9, B11, B13, B15, B17]]

	A = sp.Matrix([[A0+A1,A2+A3,A4+A5],[A6+A7,A8+A9,A10+A11],[A12+A13,A14+A15,A16+A17]])
	B = sp.Matrix([[B0+B1,B2+B3,B4+B5],[B6+B7,B8+B9,B10+B11],[B12+B13,B14+B15,B16+B17]])
	tr1 = sp.expand((A*B)[0])
	tr2 = sp.expand((A*B)[4])
	tr3 = sp.expand((A*B)[8])

	# trace = sp.trace((A*B))
	res = []
	# print trace
	for trace in [tr1,tr2,tr3]:
		traceExpanded = sp.expand(trace)
		# U_unsorted = [re.findall("([\+\-\ ]?)[^I]?([\w\*]{4,})",str(i)) for i in UExpanded]
		traceExpanded_unsorted = re.findall("([\+\-])?[ ]?([\w\*]{4,})",str(traceExpanded))
		# print traceExpanded_unsorted

		reals = []
		imags = []

		for expr,i in zip(traceExpanded_unsorted,range(len(traceExpanded_unsorted))):
			line_imags = []
			line_reals = []

			# for expr in matrix_pos:
			# print expr
			if expr[-1].split('*')[0] == "I":
				# print "*".join(expr[-1].split("*")[1:])
				line_imags.append(" ".join([" " + expr[0],"*".join(expr[-1].split("*")[1:])]))
			else:
				line_reals.append("".join(expr))


			# print expr, line_imags, line_reals
			line_imags = convertToCpp(line_imags,N=3,symbol=["A","B"])
			line_reals = convertToCpp(line_reals,N=3,symbol=["A","B"])
			imags.append(line_imags)
			reals.append(line_reals)
			# line_reals = convertToCpp(line_reals)
			# imags.append((2*i+1,line_imags))
			# reals.append((2*i,line_reals))
		# print "".join(imags) + ";"
		res.append("".join(reals) + "\n")
	print "+".join(res) + ";"



if __name__ == '__main__':
	# print_matrix_line(symbol="A",print_complex=True)
	# print_matrix_line(symbol="B",print_complex=True)
	# print_matrix(symbol="A")
	# print_matrix(symbol="B")
	get_trace_matrix()