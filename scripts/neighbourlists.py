import numpy as np

# Processors in one dimension
Nx = 4
Ny = 4
Nz = 4
Nt = 8

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
		print_break(Np, neighbours_string)
	print "\n"
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

def print_Y_neighbours():
	print "\nY-DIRECTION"
	for Np in xrange(Nx*Ny*Nz*Nt):
		neighbours_string = "|| %3d: y+1: %3d y-1: %3d" % (Np,getYPlusOne(Np),getYMinusOne(Np))
		print neighbours_string,
		print_break(Np, neighbours_string)
	print "\n"
################################################################################################################


#### Z-Direction ###############################################################################################
def getZPlusOne(Np):
	if (((Np/(Nx*Ny) + 1) % Nz) == 0):
		x_count = Np % Nx
		y_count = (Np/Nx) % Ny * Nx
		z_count = 0
		t_count = Np/(Nx*Ny*Nz) % Nt * (Nx*Ny*Nz)
		return x_count + y_count + z_count + t_count
	else:
		# return 0
		return Np + Nx*Ny

def getZMinusOne(Np):
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
	else:
		return Np + Nx*Ny*Nz

def getTMinusOne(Np):
	x_count = Np % Nx
	y_count = ((Np / Nx) % Ny) * Nx
	z_count = (Np/(Nx*Ny) % Nz) * (Nx*Ny)
	t_count = (Np/(Nx*Ny*Nz) - 1) % Nt * (Nx*Ny*Nz)
	return x_count + y_count + z_count + t_count

def print_T_neighbours():
	print "\nT-DIRECTION"
	for Np in xrange(Nt*Nz*Ny*Nx):
		neighbours_string = "|| %3d: t+1: %3d t-1: %3d" % (Np,getTPlusOne(Np),getTMinusOne(Np))
		print neighbours_string,
		print_break(Np, neighbours_string)
	print "\n"
################################################################################################################

def print_processor_lattice():
	for Np in xrange(Nz*Ny*Nx*Nt):
		print "%2d" % Np,
		if (Np+1) % (Nx*Ny*Nz) == 0:
			print "\n"+6*"="
		elif (Np+1) % (Nx*Ny) == 0:
			print "\n"
		elif (Np+1) % Nx == 0:
			print ""

def getNeighbourLists():
	n_lists = []
	for Np in xrange(Nt*Nz*Ny*Nx):
		n_lists.append([getXMinusOne(Np),getXPlusOne(Np),
						getYMinusOne(Np),getYPlusOne(Np),
						getZMinusOne(Np),getZPlusOne(Np),
						getTMinusOne(Np),getTPlusOne(Np)])
	return n_lists

def print_neighbours(n_lists,Np):
	n_list = n_lists[Np]
	print "Process rank = %d" % Np
	print "x-1 = %2d x+1 = %2d" % (n_list[0],n_list[1])
	print "y-1 = %2d y+1 = %2d" % (n_list[2],n_list[3])
	print "z-1 = %2d z+1 = %2d" % (n_list[4],n_list[5])
	print "t-1 = %2d t+1 = %2d" % (n_list[6],n_list[7])

def main():
	# ## Printing neighbour-lists
	# print_X_neighbours()
	# print_Y_neighbours()
	# print_Z_neighbours()
	# print_T_neighbours()

	neighbour_lists = getNeighbourLists()
	print_processor_lattice()

	# for i in range(Nx*Ny*Nz*Nt):
	# 	print_neighbours(neighbour_lists,i)

	# print_neighbours(neighbour_lists[0],0)

if __name__ == '__main__':
	main()