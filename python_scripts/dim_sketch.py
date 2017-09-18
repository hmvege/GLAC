dims = 8
counter = 0

def dim(n):
	if n%4==0:
		return "x"
	elif n%4==1:
		return "y"
	elif n%4==2:
		return "z"
	else:
		return "t"

for i in range(dims):
	msg = "="*8 + "new face at %s" % dim(i) + "="*8
	print msg
	for j in range(dims):
		if i!=j:
			print i,j
			counter += 1
	print len(msg) * "=" + "\n"

print counter