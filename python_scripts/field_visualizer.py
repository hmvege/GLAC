import numpy as np
import os
import itertools
from mayavi import mlab
import subprocess

def create_animation(figpath, observable, method, anim_type):
	"""
	Method for created gifs and movies from generated 3D figures.

	Args:
		figpath: path to figures.
		observable: observable we are creating a gif or a movie for.
		method: type of 3D plot.
	"""
	_ANIM_TYPES = ["gif", "mp4"]
	assert anim_type in _ANIM_TYPES, "%s is not a recognized animation type." % anim_type
	
	output_path = os.path.join(figpath, '%s_%s.%s' % (observable, method, anim_type))
	if anim_type == "gif":
		input_paths = os.path.join(figpath, '%s_%s_t*.png' % (observable, method))
		cmd = ['convert', '-delay', '1', '-loop', '0', input_paths, output_path]
	else:
		input_paths = os.path.join(figpath, '%s_%s_t%%02d.png' % (observable, method))
		cmd = ['ffmpeg', '-framerate', '4', '-i', input_paths, '-c:v', 'libx264', '-r', '30', output_path]
	subprocess.Popen(cmd, stdout=subprocess.PIPE)
	print "%s animation %s created." % (anim_type, output_path)

N = 8
NT = 16
threads = 8
observable = "energy"
flowtime = 1000

file_name = ("../output/lattice_field_density{0:<d}x{1:<d}/scalar_fields/"
			 "{3:<s}/lattice_field_density{0:<d}x{1:<d}_"
			 "b6.000000_N{0:<d}_NT{1:<d}_np{2:<d}_config0{4:<04d}.bin").format(N, NT, threads, observable, flowtime)

def global_index(i, j, k, l):
	return i + N*(j + N*(k + N*l)) # column-major

lattice_size = N**3 * NT
field = np.zeros((N, N, N, NT))

file = np.fromfile(file_name)

print "Retrieved file %s" % file_name

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

# field *= 1000

print "Populated field."

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

def plot_points3d(F, ftype, figpath, obs, vmin=None, vmax=None, cgif=True, cmovie=True):
	mlab.figure()
	for it in xrange(NT):
		mlab.clf()
		fpath = os.path.join(figpath, "%s_points3d_t%02d.%s" % (obs, it, ftype))
		# mlab.points3d(F[:,:,:,it], vmin=vmin, vmax=vmax, scale_factor=100, scale_mode="scalar")
		mlab.points3d(F[:,:,:,it], vmin=vmin, vmax=vmax)
		mlab.savefig(fpath)
		print "file created at %s" % fpath
	mlab.close()

	if cgif:
		create_animation(figpath, observable, "points3d", "gif")
	if cmovie:
		create_animation(figpath, observable, "points3d", "mp4")

# Small routine for making program sleep
import wx
import time
def animate_sleep(x):
	n_steps = int(x / 0.01)
	for i in range(n_steps):
		time.sleep(0.01)
		wx.Yield()

def plot_scalar_field(F, ftype, figpath, obs, anim=False, vmin=None, vmax=None, cgif=True, cmovie=True):
	mlab.figure()
	if anim:
		source = mlab.pipeline.scalar_field(F[:,:,:,0])
		vol = mlab.pipeline.volume(source, vmin=vmin, vmax=vmax)


	# # Changing the ctf:
	# from tvtk.util.ctf import ColorTransferFunction
	# ctf = ColorTransferFunction()
	# ctf.add_rgb_point(vmin, 0, 0, 1)  # r, g, and b are float
	# ctf.add_rgb_point(vmax, 1, 0, 0)  # r, g, and b are float
	#                                    # between 0 and 1
	# # ...
	# vol._volume_property.set_color(ctf)
	# vol._ctf = ctf
	# vol.update_ctf = True

	# # Changing the otf:
	# from tvtk.util.ctf import PiecewiseFunction
	# otf = PiecewiseFunction()
	# otf.add_point(vmin, 0)
	# otf.add_point(vmax, 1)
	# vol._otf = otf
	# vol._volume_property.set_scalar_opacity(otf)

	for it in xrange(1,NT):
		fpath = os.path.join(figpath, "%s_volume_t%02d.%s" % (obs, it, ftype))
		mlab.clf()
		
		if anim:
			animate_sleep(1)
			vol.mlab_source.set(scalars=F[:,:,:,it], vmin=vmin, vmax=vmax)
		else:
			source = mlab.pipeline.scalar_field(F[:,:,:,it])
			vol = mlab.pipeline.volume(source, vmin=vmin, vmax=vmax)

			# # Changing the ctf:
			# from tvtk.util.ctf import ColorTransferFunction
			# ctf = ColorTransferFunction()
			# ctf.add_rgb_point(vmin, 0, 0, 1)  # r, g, and b are float
			# ctf.add_rgb_point(vmax, 1, 0, 0)  # r, g, and b are float
			#                                    # between 0 and 1
			# # ...
			# vol._volume_property.set_color(ctf)
			# vol._ctf = ctf
			# vol.update_ctf = True

			# # Changing the otf:
			# from tvtk.util.ctf import PiecewiseFunction
			# otf = PiecewiseFunction()
			# otf.add_point(vmin, 0)
			# otf.add_point(vmax, 1)
			# vol._otf = otf
			# vol._volume_property.set_scalar_opacity(otf)


		mlab.savefig(fpath)
		print "file created at %s" % fpath

	mlab.close()

	if cgif:
		create_animation(figpath, observable, "volume", "gif")
	if cmovie:
		create_animation(figpath, observable, "volume", "mp4")

plot_points3d(field, file_type, animation_figure_fpath, observable, vmin=min_val*1.5, vmax=max_val*0.9)
# plot_scalar_field(field, file_type, animation_figure_fpath, observable, anim=False, vmin=min_val*1.5, vmax=max_val*0.9)

print "Done"