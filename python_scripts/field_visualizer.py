import numpy as np
import os
import itertools

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import scipy.interpolate as scp

import traits

from mayavi import mlab

N = 8
NT = 16
threads = 8
observable = "energy"
flowtime = 1000

file_name = ("../output/lattice_field_density{0:<d}x{1:<d}/scalar_fields/"
			 "{3:<s}/lattice_field_density{0:<d}x{1:<d}_"
			 "b6.000000_N{0:<d}_NT{1:<d}_np{2:<d}_config0{4:<04d}.bin").format(N, NT, threads, observable, flowtime)

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

file_type = "png"

x = np.zeros(np.prod(field[:,:,:,0].shape))
y = np.zeros(np.prod(field[:,:,:,0].shape))
z = np.zeros(np.prod(field[:,:,:,0].shape))
for i, n in enumerate(itertools.product(range(N), repeat=3)):
	x[i], y[i], z[i] = n

# scalar_array = np.zeros(np.prod(field[:,:,:,0].shape))
scalar_array = np.zeros((N, N, N))

for ix in xrange(N):
	for iy in xrange(N):
		for iz in xrange(N):
			scalar_array[ix,iy,iz] = field[ix,iy,iz,0]

max_val = np.max(field)
min_val = np.min(field)

def plot_points3d(F, ftype, figpath):
	mlab.figure()
	for it in xrange(NT):
		fpath = os.path.join(figpath, "act_dens_t%02d.%s" % (it, ftype))
		mlab.points3d(F[:,:,:,it])
		# mlab.savefig(fpath)
		mlab.show()
		# mlab.clf()
		print "file created at %s" % fpath

# plot_points3d(field, file_type, animation_figure_fpath)

# x,y,z = np.mgrid[0:N:1,0:N:1,0:N:1]

# mlab.contour3d(scalar_array, vmax=max_val, vmin=min_val, transparent=True, contours=4)
# # mlab.contour3d(x,y,z, scalar_array, vmax=max_val, vmin=min_val, transparent=True, contours=4)

# Small routine for making program sleep
import wx
import time
def animate_sleep(x):
	n_steps = int(x / 0.01)
	for i in range(n_steps):
		time.sleep(0.01)
		wx.Yield()

def plot_points3d(F, ftype, figpath):
	mlab.figure()
	# source = mlab.pipeline.scalar_field(F[:,:,:,0])
	# vol = mlab.pipeline.volume(source, vmin=min_val*2, vmax=max_val*0.8)

	for it in xrange(1,NT):
		mlab.clf()
		# animate_sleep(1)
		fpath = os.path.join(figpath, "act_dens_t%02d.%s" % (it, ftype))
		# vol.mlab_source.set(scalars=F[:,:,:,it])

		source = mlab.pipeline.scalar_field(F[:,:,:,it])
		vol = mlab.pipeline.volume(source, vmin=min_val*2, vmax=max_val*0.8)

		mlab.savefig(fpath)
		# mlab.show()
		print "file created at %s" % fpath

	mlab.close()

plot_points3d(field, file_type, animation_figure_fpath)
# mlab.show()
print "Done"


