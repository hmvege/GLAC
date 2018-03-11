import numpy as np
import os
from mayavi import mlab

file_name = ("../output/lattice_field_density_test_run/field_configurations/"
			 "lattice_field_density_test_run_energyflowLatticeDoublesField_beta6.000000_spatial8_temporal16_threads8_config00000.bin")

def global_index(i, j, k, l):
	return i + N*(j + N*(k + N*l)) # column-major

N = 8
NT = 16
lattice_size = N**3 * NT
field = np.zeros((N, N, N, NT))

file = np.fromfile(file_name)

print "Retrieved file"
# print file.shape
# print np.sum(file) / lattice_size


# with open(file_name, "rb") as f:
# 	byte = f.read(1)
# 	while byte != "":
# 		for it in xrange(NT):
# 			for iz in xrange(N):
# 				for iy in xrange(N):
# 					for ix in xrange(N):
# 						# field[ix, iy, iz, it] = file[global_index(ix, iy, iz, it)]
# 						# try:
# 						field[ix, iy, iz, it] = byte
# 						# except:
# 						# 	print byte

# print np.sum(field)

for it in xrange(NT):
	for iz in xrange(N):
		for iy in xrange(N):
			for ix in xrange(N):
				field[ix, iy, iz, it] = file[global_index(ix, iy, iz, it)]

print "Populated field."

print field[:,:,:,0].shape
# exit(1)
# Plot scatter with mayavi
# figure = mlab.figure('DensityPlot')
# pts = mlab.points3d(field[:,:,:,0], scale_mode='none', scale_factor=0.07)

# mlab.axes()

# source = mlab.pipeline.scalar_field(field[:,:,:,0])
# print source
# vol = mlab.pipeline.volume(source)
mlab.figure()
# mlab.test_contour3d()
mlab.contour3d(field[:,:,:,0])
# mlab.test_plot3d()# (field[:,:,:,0])
# mlab.points3d(field[:,:,:,0])
mlab.show()

print "Done"


