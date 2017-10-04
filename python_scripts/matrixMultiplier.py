msg = ""

N = 2

for i in range(N):
	for j in range(N):
		msg += "temp[%d] = " % (N*i+j)
		for k in range(N):
			msg += "mat[%d]*B.mat[%d]" % (N*i+k,N*k+j)
			if k != N-1: msg += " + "
		print msg + ";"
		msg = ""