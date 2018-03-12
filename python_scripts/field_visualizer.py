import numpy as np
import os

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import scipy.interpolate as scp

from mayavi import mlab

N = 4
NT = 8
threads = 2

file_name = ("../output/lattice_field_density{0:<d}x{1:<d}/field_configurations/"
			 "lattice_field_density{0:<d}x{1:<d}_energyflowLatticeDoublesField_"
			 "beta6.000000_spatial{0:<d}_temporal{1:<d}_threads{2:<d}_config01000.bin").format(N, NT, threads)

# file_name = ("../output/lattice_field_density_test_run/field_configurations/"
# 			 "lattice_field_density_test_run_topcflowLatticeDoublesField_beta6.000000_spatial8_temporal16_threads8_config01000.bin")


def global_index(i, j, k, l):
	return i + N*(j + N*(k + N*l)) # column-major

lattice_size = N**3 * NT
field = np.zeros((N, N, N, NT))

file = np.fromfile(file_name)

print "Retrieved file %s" % file_name
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
				field[ix, iy, iz, it] = - file[global_index(ix, iy, iz, it)] / 64.0

print "Populated field."

print field[:,:,:,0].shape
# exit(1)
# Plot scatter with mayavi
# figure = mlab.figure('DensityPlot')
# pts = mlab.points3d(field[:,:,:,0], scale_mode='none', scale_factor=0.07)


#### Bad matplotlib method
# fig = plt.figure()
# ax = fig.gca(projection='3d')

# for iz in xrange(N):
# 	for iy in xrange(N):
# 		for ix in xrange(N):
# 			# field[ix, iy, iz, it] = file[global_index(ix, iy, iz, it)]
# 			# try:
# 			# print field[ix, iy, iz, 0]
# 			ax.scatter(ix, iy, iz, c=field[ix, iy, iz, 0], marker="o", s=100, alpha=0.1, cmap="PRGn", edgecolor="0")
# 			# except:
# 			# 	print byte

# plt.show()

# print ax.scatter.__doc__

# print scp.interp2d.__doc__




# mlab.axes()

# source = mlab.pipeline.scalar_field(field[:,:,:,0])
# vol = mlab.pipeline.volume(source)
# mlab.figure()
# mlab.test_contour3d()
# mlab.contour3d(field[:,:,:,0])
# mlab.test_plot3d()# (field[:,:,:,0])

animation_figure_fpath = "../figures/act_density"

if not os.path.isdir(animation_figure_fpath):
	os.mkdir(animation_figure_fpath)
	print ">mkdir %s" % animation_figure_fpath

# print mlab.savefig.__doc__
# exit(1)

file_type = "png"


for it in xrange(NT):
	fpath = os.path.join(animation_figure_fpath, "act_dens_t%02d.%s" % (it, file_type))
	mlab.points3d(field[:,:,:,it])
	# mlab.savefig(fpath)
	mlab.show()
	# mlab.clf()
	print "file created at %s" % fpath


print "Done"


