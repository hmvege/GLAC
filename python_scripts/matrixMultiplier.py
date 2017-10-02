msg = ""



for i in range(3):
	for j in range(3):
		msg += "temp[%d] = " % (3*i+j)
		for k in range(3):
			msg += "mat[%d]*B.mat[%d]" % (3*i+k,3*k+j)
			if k != 2: msg += " + "
		print msg + ";"
		msg = ""