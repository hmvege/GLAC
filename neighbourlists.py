"""
	y
x	0	1	2	3
	4	5	6	7
	8	9	10	11
	12	13	14	15
"""

# Processors in one dimension
Nx = 4
Ny = 4
Nz = 4
Nt = 8

# X DIRECTION
print "Getting x neighbours for + 1 direction"
for Np in xrange(16):
	if ((Np + 1) % Nx) == 0:
		print "i = %5d. i+1=%5d" % (Np, (Np/Nx * Nx))
	else:
		print "i = %5d. i+1=%5d" % (Np, Np+1)

print "Getting x neighbours for - 1 direction"
for Np in xrange(16):
	print "i = %5d. i-1=%5d" % (Np, (Np-1)%Nx + Np/Nx * Nx)

# Y DIRECTION
print "Getting y neighbours for + 1 direction"
for i in xrange(4):
	for j in xrange(0,13,4):
		Np = i+j
		
		if ((Np/Nx+1) % Ny) == 0:
			print "i = %5d. i+1=%5d" % (Np, Np % Ny)
		else:
			print "i = %5d. i+1=%5d" % (Np, Np+Nx)

print "Getting y neighbours for - 1 direction"

for i in xrange(4):
	for j in xrange(0,13,4):
		Np = i+j
		# print "i = %5d. i-1=%5d" % (Np, (Np/Nx-1)%Nx * Ny )
		# if ( Np/(Nx*Ny) %  )
		print Np/(Nx)
		# print "i = %5d. i-1=%5d" % (Np, Np - Nx)
		# if ((Np/Nx-1) % Ny) == 0:
		# 	print "i = %5d. i+1=%5d" % (Np, Np % Ny)
		# else:
		# 	print "i = %5d. i+1=%5d" % (Np, Np+Nx)
