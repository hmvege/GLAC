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

N = Nx*Ny

def print_break(Np, neighbours_string):
	if (Np+1) % (Nx*Ny*Nz) == 0:
		print "\n\n\n"
	elif (Np+1) % (Nx*Ny) == 0:
		print "\n"+"="*(len(neighbours_string)*Nx+3)
	elif (Np+1) % Nx == 0:
		print ""


#### X-Direction ###############################################################################################
def getXPlusOne(Np):
	if ((Np + 1) % Nx) == 0:
		x_count = (Np/Nx * Nx)
		return x_count
	else:
		return Np+1

def getXMinusOne(Np):
	x_count = (Np-1) % Nx
	y_count = Np/Nx * Nx
	return x_count + y_count
	# OLD METHOD
	# return (Np-1) % Nx + Np/Nx * Nx

def print_X_neighbours():
	print "\nX-DIRECTION"
	for Np in xrange(Nx*Ny*Nz*Nt):
		# Method for +1
		neighbours_string = "|| %3d: x+1: %3d x-1: %3d" % (Np,getXPlusOne(Np),getXMinusOne(Np))
		print neighbours_string,
		
		# if (Np+1) % (Nx*Ny) == 0:
		# 	print "\n"+"="*(len(neighbours_string)*Nx+3)
		# elif (Np+1) % Nx == 0:
		# 	print ""
		print_break(Np, neighbours_string)
	print "\n"

	## OLD METHOD
	# for Np in xrange(Nx*Ny):
	# 	# Method for +1
	# 	if ((Np + 1) % Nx) == 0:
	# 		print "i = %2d, i+1 = %2d," % (Np, (Np/Nx * Nx)),
	# 	else:
	# 		print "i = %2d, i+1 = %2d," % (Np, Np+1),
	# 	# Method for -1
	# 	print "i-1 = %2d" % ((Np-1)%Nx + Np/Nx * Nx)
	# 	if (Np+1) % Nx == 0:
	# 		print ""
################################################################################################################


#### Y-Direction ###############################################################################################
def getYPlusOne(Np):
	if ((Np/Nx + 1) % Ny) == 0:
		x_count = Np % Nx
		y_count = (Np / (Ny*Nx)) * (Nx*Ny) # Method for getting base layer by doing integer division
		return y_count + x_count
	else:
		return Np + Nx

def getYMinusOne(Np):
	x_count = Np % Nx
	y_count = (Np/Nx - 1) % Ny * Nx
	z_count = Np/(Nx*Ny) * (Nx*Ny)
	return x_count + y_count + z_count
	# z_count = (Np/(Nx*Ny) - 1) % Nz * (Nx*Ny)
	# return (Np/Nx - 1) % Ny * Ny + (Np - Nx) % Ny

def print_Y_neighbours():
	print "\nY-DIRECTION"
	for Np in xrange(Nx*Ny*Nz*Nt):
		# Method for +1
		neighbours_string = "|| %3d: y+1: %3d y-1: %3d" % (Np,getYPlusOne(Np),getYMinusOne(Np))
		print neighbours_string,

		# if (Np+1) % (Nx*Ny) == 0:
		# 	print "\n"+"="*(len(neighbours_string)*Nx+3)
		# elif (Np+1) % Nx == 0:
		# 	print ""
		print_break(Np, neighbours_string)
	print "\n"

	## OLD METHOD
	# for i in xrange(4):
	# 	for j in xrange(0,13,4):
	# 		Np = i+j
	# 		# Method for +1
	# 		if ((Np/Nx+1) % Ny) == 0:
	# 			print "i = %2d, i+1 = %2d" % (Np, Np % Ny),
	# 		else:
	# 			print "i = %2d, i+1 = %2d" % (Np, Np+Nx),
	# 		# Method for -1
	# 		print "i-1 = %2d" % ((Np/Nx - 1) % Nx * Ny + (Np - Nx) % Ny)
	# 	print ""
################################################################################################################


#### Z-Direction ###############################################################################################
def getZPlusOne(Np):
	if (((Np/(Nx*Ny) + 1) % Nz) == 0):
		# return Np % Nz + (Np/Nz) % Ny * Nz
		return Np % Nx + (Np/Nx) % Ny * Nx
	else:
		return Np + Nx*Ny

def getZMinusOne(Np):
	# Np / Nx % Ny * Nx takes care of y direction counting
	# Np % Nx takes care of the x diretion counting
	# (Np/(Nx*Ny) - 1) % Nz * (Nx*Ny) takes care of the z direction counting
	# return (Np/(Nx*Ny) - 1) % Nz * (Nx*Ny) + Np / Nx % Ny * Nx + Np % Nx
	x_count = Np % Nx
	y_count = Np / Nx % Ny * Nx
	z_count = (Np/(Nx*Ny) - 1) % Nz * (Nx*Ny)
	t_count = Np/(Nx*Ny*Nz) * (Nx*Ny*Nz)
	return x_count + y_count + z_count + t_count

def print_Z_neighbours():
	print "\nZ-DIRECTION"
	for Np in xrange(Nz*Ny*Nx*Nt):
		neighbours_string = "|| %3d: z+1: %3d z-1: %3d" % (Np,getZPlusOne(Np),getZMinusOne(Np))
		print neighbours_string,
		# print "|| %2d: z+1: %2d" % (Np,getZPlusOne(Np)),

		# if (Np+1) % (Nx*Ny) == 0:
		# 	print "\n"+"="*(len(neighbours_string)*Nx+3)
		# elif (Np+1) % Nx == 0:
		# 	print ""
		print_break(Np, neighbours_string)
	print "\n"
################################################################################################################


#### T-Direction ###############################################################################################
def getTPlusOne(Np):
	if (((Np/(Nx*Ny*Nz) + 1) % Nt) == 0):
		x_count = Np % Nx
		y_count = (Np/Nx) % Ny * Nx
		z_count = (Np/(Nx*Ny)) % Nz * (Nx*Ny)
		return x_count + y_count + z_count
		# return Np % Nx + (Np/Nx) % Ny * Nx
	else:
		return Np + Nx*Ny*Nz

def getTMinusOne(Np):
	x_count = Np % Nx
	y_count = Np / Nx % Ny * Nx
	z_count = Np/(Nx*Ny) % Nz * (Nx*Ny)
	t_count = (Np/(Nx*Ny*Nz) - 1) % Nt * (Nx*Ny*Nz)
	return x_count + y_count + z_count + t_count
	# return z_count

def print_T_neighbours():
	print "\nT-DIRECTION"
	for Np in xrange(Nt*Nz*Ny*Nx):
		neighbours_string = "|| %3d: t+1: %3d t-1: %3d" % (Np,getTPlusOne(Np),getTMinusOne(Np))
		print neighbours_string,
		# print "|| %2d: z+1: %2d" % (Np,getZPlusOne(Np)),

		# if (Np+1) % (Nx*Ny*Nz) == 0:
		# 	print "\n\n\n"
		# elif (Np+1) % (Nx*Ny) == 0:
		# 	print "\n"+"="*(len(neighbours_string)*Nx+3)
		# elif (Np+1) % Nx == 0:
		# 	print ""
		print_break(Np, neighbours_string)
	print "\n"
################################################################################################################

def print_processor_lattice():
	for Np in xrange(Nz*Ny*Nx):
		print "%2d" % Np,
		if (Np+1) % (Nx*Ny) == 0:
			print "\n"
		elif (Np+1) % Nx == 0:
			print ""

# print_X_neighbours()
# print_Y_neighbours()
# print_Z_neighbours()
# print_T_neighbours()

# print_processor_lattice()