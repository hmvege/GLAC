import numpy as np
import os
import itertools
from mayavi import mlab
import subprocess
import wx
import time
import copy as cp
import re

# Small routine for making program sleep
def animate_sleep(x):
	n_steps = int(x / 0.01)
	for i in range(n_steps):
		time.sleep(0.01)
		wx.Yield()

def check_folder(folder_name, dryrun, verbose=False):
	# Checks that figures folder exist, and if not will create it
	if not os.path.isdir(folder_name):
		if dryrun or verbose:
			print "> mkdir %s" % folder_name
		if not dryrun:
			os.mkdir(folder_name)

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

	animation_figure_path = os.path.join(figpath, '%s_%s.%s' % (observable, method, anim_type))
	if anim_type == "gif":
		input_paths = os.path.join(figpath, '%s_t*.png' % method)
		cmd = ['convert', '-delay', '1', '-loop', '0', input_paths, animation_figure_path]
	else:
		input_paths = os.path.join(figpath, '%s_t%%02d.png' % method)
		cmd = ['ffmpeg', '-framerate', '4', '-i', input_paths, '-y', '-c:v', 'libx264', '-r', '30', animation_figure_path]
	subprocess.Popen(cmd, stdout=subprocess.PIPE)
	print "%s animation %s created." % (anim_type, animation_figure_path)

def global_index(i, j, k, l, N):
	return i + N*(j + N*(k + N*l)) # column-major

def plot_points3d(F, ftype, figpath, obs, NT, vmin=None, vmax=None, cgif=True, cmovie=True):
	mlab.figure()
	factor = 100.
	F *= factor
	for it in xrange(NT):
		mlab.clf()
		fpath = os.path.join(figpath, "%s_points3d_t%02d.%s" % (obs, it, ftype))

		mlab.points3d(F[:,:,:,it], vmin=vmin*factor, vmax=vmax*factor, scale_factor=1/(vmax*factor/np.max(F[:,:,:,it])), scale_mode="scalar")
		# mlab.points3d(F[:,:,:,it], vmin=vmin, vmax=vmax)
		mlab.savefig(fpath)
		print "file created at %s" % fpath
	mlab.close()

	if cgif:
		create_animation(figpath, obs, "points3d", "gif")
	if cmovie:
		create_animation(figpath, obs, "points3d", "mp4")

def plot_scalar_field(F, ftype, figpath, obs, anim=False, vmin=None, vmax=None, cgif=True, cmovie=True):
	mlab.figure()
	if anim:
		source = mlab.pipeline.scalar_field(F[:,:,:,0])
		vol = mlab.pipeline.volume(source, vmin=vmin, vmax=vmax)

	for it in xrange(1,NT):
		fpath = os.path.join(figpath, "%s_volume_t%02d.%s" % (obs, it, ftype))
		mlab.clf()
		
		if anim:
			animate_sleep(1)
			vol.mlab_source.set(scalars=F[:,:,:,it], vmin=vmin, vmax=vmax)
		else:
			source = mlab.pipeline.scalar_field(F[:,:,:,it])
			vol = mlab.pipeline.volume(source, vmin=vmin, vmax=vmax)

		mlab.savefig(fpath)
		print "file created at %s" % fpath

	mlab.close()

	if cgif:
		create_animation(figpath, obs, "volume", "gif")
	if cmovie:
		create_animation(figpath, obs, "volume", "mp4")

def plot_iso_surface(F, ftype, figpath, obs, NT, vmin=None, vmax=None, cgif=True, cmovie=True):
	mlab.figure()
	for it in xrange(NT):
		mlab.clf()
		fpath = os.path.join(figpath, "%s_iso_surface_t%02d.%s" % (obs, it, ftype))
		# mlab.points3d(F[:,:,:,it], vmin=vmin, vmax=vmax, scale_factor=100, scale_mode="scalar")
		source = mlab.pipeline.scalar_field(F[:,:,:,it])
		mlab.pipeline.iso_surface(source, vmin=vmin*1.1, vmax=vmax*0.9, contours=6)
		mlab.savefig(fpath)
		print "file created at %s" % fpath
	mlab.close()

	if cgif:
		create_animation(figpath, obs, "iso_surface", "gif")
	if cmovie:
		create_animation(figpath, obs, "iso_surface", "mp4")

class FieldAnimation:
	"""
	Class for handling the animation of lattice fields.
	"""
	def __init__(self, batch_name, N, NT, verbose=False, dryrun=False):
		"""
		Initializer for the field animator.
		"""
		self.verbose = verbose
		self.dryrun = dryrun
		self.N = N
		self.NT = NT
		self.data = {}
		self.batch_name = batch_name

		self._retrieve_field_data()

		# Creates folder to place animations in
		temp_output_folder = os.path.join("..", "figures", self.batch_name)
		check_folder(temp_output_folder, self.dryrun, verbose=self.verbose)
		self.output_folder = os.path.join(temp_output_folder, "field_animations")
		check_folder(self.output_folder, self.dryrun, verbose=self.verbose)
		for key in self.data:
			obs_folder_paths = os.path.join(self.output_folder, key)
			check_folder(obs_folder_paths, self.dryrun, verbose=self.verbose)
			self.data[key]["animation_figure_path"] = obs_folder_paths

	def _retrieve_field_data(self):
		"""Function for retrieving fields."""

		# Load the different data of the fields
		scalar_field_folder = os.path.join("..", 
			"output", self.batch_name, "scalar_fields")
		
		# Goes through the observables in the scalar fields folder
		for obs in os.listdir(scalar_field_folder):
			obs_folder = os.path.join(scalar_field_folder, obs)

			# Ensures we have flow time files in the observable folder
			if len(os.listdir(obs_folder)) == 0:
				continue 
			
			self.data[obs] = cp.deepcopy(self._get_flow_files(obs_folder))
			
			# Factors in missing -1.0/64.0
			if obs == "energy":
				for t in self.data[obs].keys():
					if not isinstance(t, int):
						continue
					self.data[obs][t] *= (-1.0/64.0)

			if self.verbose:
				print "Data retrieved for observable %s." % obs

	def _get_flow_files(self, obs_folder):
		"""
		Function for populating a data dicitonary for a given observable folder.
		"""
		data_dict = {}

		# Goes through flow observables in observable folder
		for flow_obs_file in os.listdir(obs_folder):
			flow_obs_file_path = os.path.join(obs_folder, flow_obs_file)
			
			# Gets data array
			raw_data = np.fromfile(flow_obs_file_path)

			if self.verbose:
				print "Retrieved file %s" % flow_obs_file_path

			# Gets the flow time
			flow_time = self._get_cfg_num_from_file(flow_obs_file)

			# Converts data array to data field, and sets it in the dictionary
			data_dict[flow_time] = self._convert_to_field(raw_data) 

		return data_dict

	def _convert_to_field(self, data_array):
		"""Converts data array to a field."""
		_field = np.zeros((self.N, self.N, self.N, self.NT))
		for it in xrange(self.NT):
			for iz in xrange(self.N):
				for iy in xrange(self.N):
					for ix in xrange(self.N):
						_field[ix, iy, iz, it] = data_array[self._global_index(ix, iy, iz, it)]

		return _field

	def _global_index(self, i, j, k, l):
		"""Column-major contigious memory locator."""
		return int(i + self.N*(j + self.N*(k + self.N*l)))

	def _get_cfg_num_from_file(self, fname):
		"""
		Function for getting the flow time from a file name by looking 
		for 'config'.

		Args:
			fname: string file name to look for config number in.
		
		Returns:
			integer config number.
		"""
		_flow_time_list = re.findall("config(\d+)", fname)
		assert len(_flow_time_list) == 1, "error in file name: %s" % fname
		return int(_flow_time_list[0])

	def _get_output_animation_folder(self, output_folder):
		if output_folder != None:
			return output_folder
		else:
			return self.data[observable]["animation_figure_path"]


	def animate(self, observable, time_type, time_slice, plot_type, **kwargs):
		"""
		Method for animating in flow time or euclidean time.

		Args:
			observable: string name of observable we are looking at.
			time_type: "euclidean" or "flow". Will plot evolution in provided time type.
			time_slice: integer in euclidean time we are looking at.
			plot_type: type of plot we are having.
			**kwargs passed to the animation function.
		"""

		# Sets up all of the available flow time in a sorted list
		flow_times = sorted([key for key in self.data[observable].keys() if isinstance(key, int)])

		# Sets up folders
		time_type_folder_path = os.path.join(self.data[observable]["animation_figure_path"], time_type)
		check_folder(time_type_folder_path, self.dryrun, verbose=self.verbose)

		if time_type == "flow":
			# For plotting evolution in flow time

			field_data = []
			
			if time_slice >= self.data[observable][0].shape[-1]:
				raise IndexError(("Out of bounds for plotting flow at Euclidean time point %d"
					" in data with points %d" % (time_slice, self.data[observable][0].shape[-1])))

			# Sets up data to be plotted			
			for t in flow_times:
				field_data.append(self.data[observable][t][:,:,:,time_slice])
			field_data = np.asarray(field_data)
			field_data = np.rollaxis(field_data, 0, 4)

			n_time_points = len(flow_times)

			time_slice_folder_path = os.path.join(time_type_folder_path, "t_eucl" + str(time_slice))
		elif time_type == "euclidean":
			# For plotting evolution in euclidean time

			if time_slice not in self.data[observable].keys():
				raise IndexError(("Out of bounds for plotting Euclidean time "
					"evolution at flow time %d with available points as %s" %
					 (time_slice, ", ".join(self.data[observable][t].keys()))))

			field_data = cp.deepcopy(self.data[observable][time_slice])
			n_time_points = self.data[observable][time_slice].shape[-1]

			time_slice_folder_path = os.path.join(time_type_folder_path, "t_flow" + str(time_slice))
		else:
			raise KeyError("Cannot plot in %s." % time_type)

		# Checks that the sub time slice folder exists
		check_folder(time_slice_folder_path, self.dryrun, verbose=self.verbose)

		# field_data = np.log(field_data)

		max_val = np.max(field_data)
		min_val = np.min(field_data)

		if plot_type == "iso_surface":
			self._plot_iso_surface(field_data, n_time_points, observable, 
				vmin=min_val, vmax=max_val, output_folder=time_slice_folder_path, **kwargs)
		elif plot_type == "volume":
			self._plot_scalar_field(field_data, n_time_points, observable, 
				vmin=min_val, vmax=max_val, output_folder=time_slice_folder_path, **kwargs)
		elif plot_type == "points3d":
			self._plot_points3d(field_data, n_time_points, observable, 
				vmin=min_val, vmax=max_val, output_folder=time_slice_folder_path, **kwargs)
		else:
			raise KeyError("Plot type %s not recognized" % plot_type)


	def _plot_iso_surface(self, F, n_time_points, observable, file_type="png", vmin=None, vmax=None, output_folder=None, cgif=True, cmovie=True):
		"""
		Function for plotting iso surfaces.

		Args:
			F: field array of size (N,N,N,NT) to plot.
			n_time_points: integer of time points NT to animate over.
			observable: str of observable we are plotting.
			file_type: string of file extension type. Default is 'png'.
			vmin: float lower cutoff value of the field. Default is None.
			vmax: float upper cutoff value of the field. Default is None.
			output_folder: location of where to place output files.
			cgif: bool if we are to create a gif.
			cmovie: bool if we are to create a movie.
		"""

		animation_figure_path = self._get_output_animation_folder(output_folder)

		mlab.figure()
		for it in xrange(n_time_points):
			mlab.clf()
			fpath = os.path.join(animation_figure_path, "iso_surface_t%02d.%s" % (it, file_type))
			# mlab.points3d(F[:,:,:,it], vmin=vmin, vmax=vmax, scale_factor=100, scale_mode="scalar")
			source = mlab.pipeline.scalar_field(F[:,:,:,it])
			mlab.pipeline.iso_surface(source, vmin=vmin, vmax=vmax, contours=4)
			mlab.savefig(fpath)
			print "file created at %s" % fpath
		mlab.close()

		if cgif:
			create_animation(animation_figure_path, observable, "iso_surface", "gif")
		if cmovie:
			create_animation(animation_figure_path, observable, "iso_surface", "mp4")

	def _plot_scalar_field(self, F, n_time_points, observable, file_type="png", vmin=None, vmax=None, output_folder=None, cgif=True, cmovie=True, live_animation=False):
		"""
		Function for plotting a scalar field.

		Args:
			F: field array of size (N,N,N,NT) to plot.
			n_time_points: integer of time points NT to animate over.
			observable: str of observable we are plotting.
			file_type: string of file extension type. Default is 'png'.
			vmin: float lower cutoff value of the field. Default is None.
			vmax: float upper cutoff value of the field. Default is None.
			output_folder: location of where to place output files.
			cgif: bool if we are to create a gif.
			cmovie: bool if we are to create a movie.
			live_animation: if we are to animate on the fly.
		"""

		animation_figure_path = self._get_output_animation_folder(output_folder)

		mlab.figure()
		if live_animation:
			source = mlab.pipeline.scalar_field(F[:,:,:,0])
			vol = mlab.pipeline.volume(source, vmin=vmin, vmax=vmax)
			start_point = 1
		else:
			start_point = 0

		for it in xrange(start_point, n_time_points):
			fpath = os.path.join(animation_figure_path, "volume_t%02d.%s" % (it, file_type))
			mlab.clf()

			if live_animation:
				animate_sleep(1)
				vol.mlab_source.set(scalars=F[:,:,:,it], vmin=vmin, vmax=vmax)
			else:
				source = mlab.pipeline.scalar_field(F[:,:,:,it])
				vol = mlab.pipeline.volume(source, vmin=vmin, vmax=vmax)

			mlab.savefig(fpath)
			print "file created at %s" % fpath
		mlab.close()

		if cgif:
			create_animation(animation_figure_path, observable, "volume", "gif")
		if cmovie:
			create_animation(animation_figure_path, observable, "volume", "mp4")

	def _plot_points3d(self, F, n_time_points, observable, file_type="png", vmin=None, vmax=None, output_folder=None, cgif=True, cmovie=True):
		"""
		Function for plotting points3d of the field.

		Args:
			F: field array of size (N,N,N,NT) to plot.
			n_time_points: integer of time points NT to animate over.
			observable: str of observable we are plotting.
			file_type: string of file extension type. Default is 'png'.
			vmin: float lower cutoff value of the field. Default is None.
			vmax: float upper cutoff value of the field. Default is None.
			output_folder: location of where to place output files.
			cgif: bool if we are to create a gif.
			cmovie: bool if we are to create a movie.
		"""

		animation_figure_path = self._get_output_animation_folder(output_folder)

		mlab.figure()
		factor = 1
		const = 0.5 # Lower limit on spheres
		F *= factor
		for it in xrange(n_time_points):
			mlab.clf()
			fpath = os.path.join(animation_figure_path, "points3d_t%02d.%s" % (it, file_type))

			# scale_factor = np.min(F[:,:,:,it]) / (factor*vmin) + (vmax - vmin) / 0.75 * factor
			# scale_factor=np.max(F[:,:,:,it])/(vmax*factor)
			scale_factor = (np.min(F[:,:,:,it]) - vmin) / (vmax - vmin) * (1 - const)  + const

			mlab.points3d(F[:,:,:,it], vmin=vmin*factor, vmax=vmax*factor, scale_factor=scale_factor, scale_mode="scalar")
			# mlab.points3d(F[:,:,:,it], vmin=vmin*factor, vmax=vmax*factor, scale_mode="scalar")
			# mlab.points3d(F[:,:,:,it], vmin=vmin, vmax=vmax)
			mlab.savefig(fpath)
			print "file created at %s" % fpath
		mlab.close()

		if cgif:
			create_animation(animation_figure_path, observable, "points3d", "gif")
		if cmovie:
			create_animation(animation_figure_path, observable, "points3d", "mp4")

def main():
	N = 8
	NT = 16
	threads = 8
	observable = "energy"
	flowtime = 1000
	verbose = True
	dryrun = False

	B60FieldAnimation = FieldAnimation("lattice_field_density8x16", N, NT, verbose=verbose, dryrun=dryrun)
	# B60FieldAnimation.animate("energy", "euclidean", 0, "iso_surface")
	# B60FieldAnimation.animate("energy", "euclidean", 400, "iso_surface")
	# B60FieldAnimation.animate("energy", "euclidean", 0, "volume")
	# B60FieldAnimation.animate("energy", "euclidean", 400, "volume")
	B60FieldAnimation.animate("energy", "euclidean", 0, "points3d")
	B60FieldAnimation.animate("energy", "euclidean", 400, "points3d")
	
	# B60FieldAnimation.animate("energy", "flow", 0, "iso_surface")
	B60FieldAnimation.animate("energy", "flow", 15, "points3d")
	# B60FieldAnimation.animate("energy", "flow", 15, "volume")

	exit(1)

	file_name = ("../output/lattice_field_density{0:<d}x{1:<d}/scalar_fields/"
				 "{3:<s}/lattice_field_density{0:<d}x{1:<d}_"
				 "b6.000000_N{0:<d}_NT{1:<d}_np{2:<d}_config0{4:<04d}.bin").format(N, NT, threads, observable, flowtime)

	lattice_size = N**3 * NT
	field = np.zeros((N, N, N, NT))

	file = np.fromfile(file_name)

	print "Retrieved file %s" % file_name

	for it in xrange(NT):
		for iz in xrange(N):
			for iy in xrange(N):
				for ix in xrange(N):
					field[ix, iy, iz, it] = - file[global_index(ix, iy, iz, it, N)] / 64.0

	print "Populated field."

	animation_figure_fpath = "../figures/act_density"

	if not os.path.isdir(animation_figure_fpath):
		os.mkdir(animation_figure_fpath)
		print ">mkdir %s" % animation_figure_fpath

	file_type = "png"

	max_val = np.max(field)
	min_val = np.min(field)

	# plot_points3d(field, file_type, animation_figure_fpath, observable, NT, vmin=min_val*1.5, vmax=max_val*0.9)

	# plot_scalar_field(field, file_type, animation_figure_fpath, observable, anim=False, vmin=min_val*1.5, vmax=max_val*0.9)

	plot_iso_surface(field, file_type, animation_figure_fpath, observable, NT, vmin=min_val, vmax=max_val)

	print "Done"


if __name__ == '__main__':
	main()