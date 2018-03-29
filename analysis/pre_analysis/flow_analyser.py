from tools.folderreadingtools import write_data_to_file, check_folder, write_raw_analysis_to_file, DataReader
from statistics.jackknife import Jackknife
from statistics.bootstrap import Bootstrap
from statistics.autocorrelation import Autocorrelation
import statistics.parallel_tools as ptools
import os
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import copy
import types

# from tqdm import trange

__all__ = ["FlowAnalyser", "AnalysePlaquette", "AnalyseEnergy",
	"AnalyseTopologicalCharge", "AnalyseTopologicalSusceptibility", 
	"AnalyseTopChargeInEuclTime", "AnalyseTopSusInEuclTime", "AnalyseQQuartic",
	"AnalyseQtQZero", "AnalyseTopChargeSplitEuclideanTime",
	"AnalysisTopSusSplitEuclideanTime", "AnalysisTopChargeSplitMCTime", 
	"AnalysisTopSusSplitMCTime", "AnalyseQtQ0Euclidean"]

class FlowAnalyser(object):
	observable_name = "Missing_Observable_Name"
	observable_name_compact = "missing_obs_name"
	x_label = "Missing x-label"
	y_label = "Missing y-label"
	mark_interval = 5
	error_mark_interval = 5
	autocorrelations_limits = 1
	folder_locations = ".."
	figures_folder = "figures"
	fname_addon = ""
	lattice_sizes = {6.0: 24**3*48, 6.1: 28**3*56, 6.2: 32**3*64, 6.45: 48**3*96}
	hbarc = 0.19732697 #eV micro m

	# Function derivative to be used in the autocorrelation class
	function_derivative = None

	# Resolution in figures created
	dpi = 350

	# Hline and Vline plots. Internal variables, accesible through children.
	plot_hline_at = None
	plot_vline_at = None

	def __init__(self, data, dryrun=False, parallel=False, numprocs=4,
		verbose=False, figures_folder=False):
		"""
		Parent class for analyzing flowed observables.

		Args:
			data: DataReader([observable_name]), an DataReader object called 
				with the compact observable name. Options: "plaq", "energy", 
				"topc".
			dryrun: optional dryrun mode. Default is False.
			parallel: optinal argument if we are to run analysis in parallel.
				Default is False.
			numprocs: optional argument for the number of processors to use. 
				Default is 4.
			verbose: optional argument for a more verbose run. Default is 
				False.
			figures_folder: optional argument for where to place the figures 
				created. Default is "../figures".

		Returns:
			Object for analyzing flow.
		"""

		# Retrieves data from data
		self.batch_name = data["batch_name"]
		self.batch_data_folder = data["batch_data_folder"]
		self.x = data["t"]
		self.y = data["obs"]
		self.flow_epsilon = data["FlowEpsilon"]

		# Sets lattice parameters
		self.beta = data["beta"]
		self.a = self.getLatticeSpacing(self.beta)
		self.r0 = 0.5 # Sommer Parameters

		# Sets up the function derivative parameters as an empty dictionary
		self.function_derivative_parameters = {}

		# Initializes up global constants
		self.N_bs = None
		self.dryrun = dryrun
		self.verbose = verbose
		if figures_folder != False:
			self.figures_folder = figures_folder
		else:
			# Adds current folder location to figures
			self.figures_folder = os.path.join(self.folder_locations, self.figures_folder)

		# Parallel variables
		self.parallel = parallel
		self.numprocs = numprocs

		# Checks that a figures folder exists
		check_folder(self.figures_folder, self.dryrun, verbose=self.verbose)

		# Check that a data run folder exist, so one data anlysis performed on different data sets do not mix
		self.data_batch_folder_path = os.path.join(self.figures_folder,os.path.split(self.batch_data_folder)[-1])
		check_folder(self.data_batch_folder_path, self.dryrun, verbose=self.verbose)

		# Checks that a batch folder exists
		self.batch_name_folder_path = os.path.join(self.data_batch_folder_path, self.batch_name)
		check_folder(self.batch_name_folder_path, self.dryrun, verbose=self.verbose)

		# Checks that observable output folder exist, and if not will create it
		self.observable_output_folder_path = os.path.join(self.batch_name_folder_path, self.observable_name_compact)
		check_folder(self.observable_output_folder_path, self.dryrun, verbose=self.verbose)

		# Sets up the post analysis folder, but do not create it till its needed
		self.post_analysis_folder_base = os.path.join(self.folder_locations, self.batch_data_folder, self.batch_name, "post_analysis_data")
		check_folder(self.post_analysis_folder_base, self.dryrun, verbose=self.verbose)

		# Checks that {post_analysis_folder}/{observable_name} exists
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_base, self.observable_name_compact)
		check_folder(self.post_analysis_folder, self.dryrun, verbose=self.verbose)

		# Makes a backup, for later use
		self.post_analysis_folder_old = self.post_analysis_folder

		# Checks if we already have scaled the x values or not
		if np.all(np.abs(np.diff(self.x) - self.flow_epsilon) > 1e-14):
			self.x = self.x * self.flow_epsilon
			self.pre_scale = False
		else:
			self.pre_scale = True

		# Max plotting window variables
		self.y_limits = [None,None]

		# Default type of observables, one per configuration per flow
		self.N_configurations, self.NFlows = self.y.shape[:2]

		# Non-bootstrapped data
		self.unanalyzed_y = np.zeros(self.NFlows)
		self.unanalyzed_y_std = np.zeros(self.NFlows)
		self.unanalyzed_y_data = np.zeros((self.NFlows, self.N_configurations))

		# Bootstrap data
		self.bootstrap_performed = False
		self.bs_y = np.zeros(self.NFlows)
		self.bs_y_std = np.zeros(self.NFlows)

		# Jackknifed data
		self.jackknife_performed = False
		self.jk_y = np.zeros(self.NFlows)
		self.jk_y_std = np.zeros(self.NFlows)

		# Autocorrelation data
		self.autocorrelation_performed = False
		self.autocorrelations = np.zeros((self.NFlows, self.N_configurations/2))
		self.autocorrelations_errors = np.zeros((self.NFlows, self.N_configurations/2))
		self.integrated_autocorrelation_time = np.ones(self.NFlows)
		self.integrated_autocorrelation_time_error = np.zeros(self.NFlows)
		self.autocorrelation_error_correction = np.ones(self.NFlows)

	def __check_ac(self, fname):
		"""If autocorrelation has been performed it will add "_noErrorCorrection" to the filename"""
		head, ext = os.path.splitext(fname)
		fname_addon = "_noErrorCorrection"
		if not self.autocorrelation_performed:
			return head + "_noErrorCorrection" + ext
		else:
			return head + ext
			
	@staticmethod
	def getLatticeSpacing(beta):
		"""
		Static method for calculating the lattice spacing a from
		the beta-value. Based on paper by M. Guagnelli et al(1998).

		Args:
			beta: beta value

		Returns:
			a: lattice spacing
		"""
		if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
		r = 0.5
		bval = (beta - 6)
		a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
		return a

	def save_post_analysis_data(self):
		"""Saves post analysis data to a file."""
		if self.bootstrap_performed:
			write_data_to_file(self)

	def save_raw_analysis_data(self, data, analysis_type):
		"""Saves the raw analysis data to post analysis folder as binary."""
		write_raw_analysis_to_file(data, analysis_type,
			self.observable_name_compact, self.post_analysis_folder,
			dryrun=self.dryrun, verbose=self.verbose)

	def boot(self, N_bs, F=None, F_error=None, store_raw_bs_values=True):
		"""
		Bootstrap caller for the flow analysis.

		Args:
			N_bs: number of bootstraps to perform.
			F: optional argument for function that will modify data. Default
				is None.
			F_error: optional argument for F function for propagating error.
				Default is None.
			store_raw_bs_values: optional argument for storing raw bootstrap
				datasets, default is True.
		"""

		# Sets default parameters
		if F == None:
			F = ptools._default_return
		if F_error == None:
			F_error = ptools._default_error_return

		# Stores number of bootstraps
		self.N_bs = N_bs

		# Sets up raw bootstrap and unanalyzed data array
		self.bs_y_data = np.zeros((self.NFlows, self.N_bs))
		self.unanalyzed_y_data = np.zeros((self.NFlows, self.N_configurations))
	
		# Generates random lists to use in resampling
		index_lists = np.random.randint(self.N_configurations, size=(self.N_bs, self.N_configurations))

		if self.parallel:
			# Sets up jobs for parallel processing
			input_values = zip(	[self.y[:,i] for i in xrange(self.NFlows)],
								[N_bs for i in xrange(self.NFlows)],
								[index_lists for i in xrange(self.NFlows)])
			
			# Initializes multiprocessing
			pool = multiprocessing.Pool(processes=self.numprocs)								

			# Runs parallel processes. Can this be done more efficiently?
			results = pool.map(ptools._bootstrap_parallel_core, input_values)

			# Garbage collection for multiprocessing instance
			pool.close()

			# Populating bootstrap data
			for i in xrange(self.NFlows):
				self.bs_y[i] = results[i][0]
				self.bs_y_std[i] = results[i][1]
				self.unanalyzed_y[i] = results[i][2]
				self.unanalyzed_y_std[i] = results[i][3]

				# Stores last data for plotting in histogram later and post analysis
				self.bs_y_data[i] = results[i][4]
				self.unanalyzed_y_data[i] = results[i][5]
		else:
			if self.verbose:
				print "Not running parallel bootstrap for %s!" % self.observable_name

			# Non-parallel method for calculating bootstrap
			for i in xrange(self.NFlows):
				bs = Bootstrap(self.y[:,i], N_bs, index_lists=index_lists)
				self.bs_y[i] = bs.bs_avg
				self.bs_y_std[i] = bs.bs_std
				self.unanalyzed_y[i] = bs.avg_original
				self.unanalyzed_y_std[i] = bs.std_original

				# Stores last data for plotting in histogram later
				self.bs_y_data[i] = bs.bs_data
				self.unanalyzed_y_data[i] = bs.data_original

		# Runs bs and unanalyzed data through the F and F_error
		self.bs_y_std = F_error(self.bs_y, self.bs_y_std)
		self.bs_y = F(self.bs_y)

		self.unanalyzed_y_std = F_error(self.unanalyzed_y, self.unanalyzed_y_std)
		self.unanalyzed_y = F(self.unanalyzed_y)

		# Runs bs data through function F
		self.bs_y_data = F(self.bs_y_data)
		self.unanalyzed_y_data = F(self.unanalyzed_y_data)

		if store_raw_bs_values:
			self.save_raw_analysis_data(self.bs_y_data, "bootstrap")

		# Sets performed flag to true
		self.bootstrap_performed = True

	def jackknife(self, F=None, F_error=None, store_raw_jk_values=True):
		"""
		When called, performs either a parallel or non-parallel jackknife, using the Jackknife class.
		
		Args:
			F: optional argument for function that will modify data. Default is None.
			F_error: optional argument for F function for propagating error. Default is None.
			store_raw_jk_values: optional argument for storing raw jackknifed datasets, default is True.

		Raises:
			AssertError will be made if F and F_error are not None or of types.FunctionType.
		"""

		# Sets default parameters
		if F == None:
			F = ptools._default_return
		if F_error == None:
			F_error = ptools._default_error_return

		# Sets up raw jackknife
		self.jk_y_data = np.zeros((self.NFlows, self.N_configurations))

		if self.parallel:
			# Sets up jobs for parallel processing
			input_values = [self.y[:,i] for i in xrange(self.NFlows)]

			# Initializes multiprocessing
			pool = multiprocessing.Pool(processes=self.numprocs)								

			# Runs parallel processes. Can this be done more efficiently?
			results = pool.map(ptools._jackknife_parallel_core, input_values)

			# Closes multiprocessing instance for garbage collection
			pool.close()

			# Populating jackknife results
			for i in xrange(self.NFlows):
				self.jk_y[i] = results[i][0]
				self.jk_y_std[i] = results[i][1]
				self.jk_y_data[i] = results[i][2]
		else:
			if self.verbose:
				print "Not running parallel jackknife for %s!" % self.observable_name

			# Non-parallel method for calculating jackknife
			for i in xrange(self.NFlows):
				jk = Jackknife(self.y[:,i])
				self.jk_y[i] = jk.jk_avg
				self.jk_y_std[i] = jk.jk_std
				self.jk_y_data[i] = jk.jk_data

		# Runs data through the F and F_error
		self.jk_y_std = F_error(self.jk_y, self.jk_y_std)
		self.jk_y = F(self.jk_y)
		self.jk_y_data = F(self.jk_y_data)

		if store_raw_jk_values:
			self.save_raw_analysis_data(self.jk_y_data, "jackknife")

		# Sets performed flag to true
		self.jackknife_performed = True

	def autocorrelation(self, store_raw_ac_error_correction=True, method="wolff"):
		"""
		Function for running the autocorrelation routine.

		Args:
			store_raw_ac_error_correction: optional argument for storing the 
				autocorrelation error correction to file.
			method: type of autocorrelation to be performed. Choose from: 
				"wolff"(default), "luscher".

		Raises:
			KeyError: if method is not a valid one.
		"""

		available_ac_methods = ["wolff", "luscher"]
		if method not in available_ac_methods:
			raise KeyError("%s not a receognized method. Choose from: %s." % (method,", ".join(available_ac_methods)))

		print "Running autocorrelation with %s ac-method" % method

		# Gets autocorrelation
		if self.parallel:
			# Sets up parallel job
			pool = multiprocessing.Pool(processes=self.numprocs)
			
			if method == "wolff":
				# Sets up jobs for parallel processing
				input_values = zip(	[self.y[:,i] for i in xrange(self.NFlows)],
									[self.function_derivative for i in xrange(self.NFlows)],
									[self.function_derivative_parameters for i in xrange(self.NFlows)])

				# Initiates parallel jobs
				results = pool.map(ptools._autocorrelation_propagated_parallel_core, input_values)
			else:
				# Sets up jobs for parallel processing
				input_values = zip([self.y[:,i] for i in xrange(self.NFlows)])
				
				# Initiates parallel jobs
				results = pool.map(ptools._autocorrelation_parallel_core, input_values)

			# Closes multiprocessing instance for garbage collection
			pool.close()

			# Populating autocorrelation results
			for i in xrange(self.NFlows):
				self.autocorrelations[i] = results[i][0]
				self.autocorrelations_errors[i] = results[i][1]
				self.integrated_autocorrelation_time[i] = results[i][2]
				self.integrated_autocorrelation_time_error[i] = results[i][3]
				self.autocorrelation_error_correction[i] = np.sqrt(2*self.integrated_autocorrelation_time[i])

		else:
			if self.verbose:
				print "Not running parallel autocorrelation for %s!" % self.observable_name

			# Non-parallel method for calculating autocorrelation
			for i in xrange(self.NFlows):
				if method == "wolff":
					ac = PropagatedAutocorrelation(self.y[:,i], 
						function_derivative=self.function_derivative, 
						func_params=self.function_derivative_parameters)
				else:
					ac = Autocorrelation(self.y[:,i])
				self.autocorrelations[i] = ac.R
				self.autocorrelations_errors[i] = ac.R_error
				self.integrated_autocorrelation_time[i] = ac.integrated_autocorrelation_time()
				self.integrated_autocorrelation_time_error[i] = ac.integrated_autocorrelation_time_error()
				self.autocorrelation_error_correction[i] = np.sqrt(2*self.integrated_autocorrelation_time[i])

		# Stores the ac error correction
		if store_raw_ac_error_correction:
			self.save_raw_analysis_data(self.autocorrelation_error_correction,
				"autocorrelation")

		# Sets performed flag to true
		self.autocorrelation_performed = True

	def plot_jackknife(self, x=None, correction_function=lambda x: x):
		"""
		Function for plotting the jackknifed data.

		Args:
			x: optional values to plot along the x-axis. Default is plotting
				self.a*np.sqrt(8*t).
			correction_function: function to correct y-axis values with.
				Default is to leave them unmodified.

		Raises:
			ValueError: if jackknife has not been performed yet.
		"""

		# Checks that jacknifing has been performed.
		if not self.jackknife_performed:
			raise ValueError("Jackknifing has not been performed yet.")

		# Sets up the x axis array to be plotted
		if isinstance(x,types.NoneType):
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x)

		# Copies data over to arrays to be plotted
		y = self.jk_y
		y_std = self.jk_y_std*self.autocorrelation_error_correction

		# Sets up the title and filename strings
		title_string = r"Jacknife of %s" % self.observable_name
		fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_jackknife_beta{1:<s}{2:<s}.png".format(self.observable_name_compact, str(self.beta).replace('.','_'), self.fname_addon))

		# Plots the jackknifed data
		self.__plot_error_core(x, correction_function(y), correction_function(y_std), title_string, fname_path)

	def plot_boot(self, x=None, correction_function=lambda x: x, _plot_bs=True):
		"""
		Function for plotting the bootstrapped data.

		Args:
			x: optional values to plot along the x-axis. Default is plotting
				self.a*np.sqrt(8*t).
			correction_function: function to correct y-axis values with.
				Default is to leave them unmodified.
			_plot_bs: internval argument for plotting the bootstrapped and
				non-bootstrapped values. Default is True.

		Raises:
			ValueError: if bootstrap has not been performed yet.
		"""

		# Checks that the bootstrap has been performed.
		if not self.bootstrap_performed and _plot_bs:
			raise ValueError("Bootstrap has not been performed yet.")

		# Retrieves relevant data and sets up the arrays to be plotted
		if isinstance(x,types.NoneType):
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x)

		# Determines if we are to plot bootstrap or original and retrieves data
		if _plot_bs:
			y = self.bs_y
			y_std = self.bs_y_std*self.autocorrelation_error_correction
		else:
			y = self.unanalyzed_y
			y_std = self.unanalyzed_y_std*self.autocorrelation_error_correction

		# Sets up the title and filename strings
		if _plot_bs:
			title_string = r"%s, $N_{bootstraps}=%d$" % (self.observable_name, self.N_bs)
			fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_bootstrap_Nbs{2:<d}_beta{1:<s}{3:<s}.png".format(self.observable_name_compact, str(self.beta).replace('.','_'), self.N_bs, self.fname_addon))
		else:
			title_string = r"%s" % self.observable_name
			fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_original_beta{1:<s}{2:<s}.png".format(self.observable_name_compact, str(self.beta).replace('.','_'), self.fname_addon))

		# Plots either bootstrapped or regular stuff
		self.__plot_error_core(x, correction_function(y), correction_function(y_std), title_string, fname_path)

	def plot_original(self, x=None, correction_function=lambda x: x):
		"""
		Plots the default analysis, mean and std of the observable.

		Args:
			x: optional values to plot along the x-axis. Default is plotting
				self.a*np.sqrt(8*t).
			correction_function: function to correct y-axis values with.
				Default is to leave them unmodified.

		"""
		self.plot_boot(x=x, correction_function=correction_function, _plot_bs=False)

	def __plot_error_core(self, x, y, y_std, title_string, fname):
		"""
		Interternal plotter function. Forms basis for unanalyzed, jackknifed 
		and bootstrapped plotters.

		Args:
			x: values to plot along the x-axis.
			y: values to plot along the y-axis.

		"""
		# Plots the jackknifed data
		fig = plt.figure()
		ax = fig.add_subplot(111)

		# Plots the error bar
		ax.errorbar(x, y, yerr=y_std, fmt=".", color="0", ecolor="r",
			label=self.observable_name, markevery=self.mark_interval,
			errorevery=self.error_mark_interval)

		# Plots hline/vline at specified position.
		# Needed for plotting at e.g. t0. 
		if self.plot_hline_at != None:
			ax.axhline(self.plot_hline_at, linestyle="--", color="0", alpha=0.3)
		if self.plot_vline_at != None:
			ax.axvline(self.plot_vline_at, linestyle="--", color="0", alpha=0.3)

		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)
		ax.set_ylim(self.y_limits)
		ax.grid(True)
		ax.set_title(title_string)
		if not self.dryrun: 
			fig.savefig(self.__check_ac(fname), dpi=self.dpi)
		if self.verbose:
			print "Figure created in %s" % fname

		# plt.show()
		plt.close(fig)

	def plot_autocorrelation(self, flow_time_index):
		"""
		Plots the autocorrelation at a given flow time.

		Args:
			flow_time_index: integer of the flow time index.

		Raises:
			ValueError: if no autocorrelation has been performed yet.
			AssertionError: if the provided flow_time_index is out of bounds.
		"""

		# Checks that autocorrelations has been performed.
		if not self.autocorrelation_performed:
			raise ValueError("Autocorrelation has not been performed yet.")

		# sets up the autocorrelation
		N_autocorr = self.N_configurations / 2

		# Converts flow_time_index if it is minus 1
		if flow_time_index == -1:
			flow_time_index = self.NFlows - 1

		# Ensures flow time is within bounds.
		assert flow_time_index < self.NFlows, "Flow time %d is out of bounds." % flow_time_index

		# Finds the maximum value at each MC time and sets up the y array
		x = range(N_autocorr)
		y = self.autocorrelations[flow_time_index,:]		
		y_std = self.autocorrelations_errors[flow_time_index,:]

		# Sets up the title and filename strings
		title_string = r"Autocorrelation of %s, $t_{flow}=%.2f$" % (self.observable_name, flow_time_index*self.flow_epsilon)
		fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_autocorrelation_flowt{1:<d}_beta{2:<s}{3:<s}.png".format(self.observable_name_compact, flow_time_index, str(self.beta).replace('.','_'), self.fname_addon))

		# Plots the autocorrelations
		fig = plt.figure()
		ax = fig.add_subplot(111)
		# ax.plot(x,y,color="0",label=self.observable_name)
		ax.errorbar(x, y, yerr=y_std, color="0", ecolor="r")#,label=self.observable_name)
		ax.set_ylim(-self.autocorrelations_limits, self.autocorrelations_limits)
		ax.set_xlim(0, N_autocorr)
		ax.set_xlabel(r"Lag $h$")
		ax.set_ylabel(r"$R = \frac{C_h}{C_0}$")
		ax.set_title(title_string)
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.2))
		ax.grid(True)
		# ax.legend()
		if not self.dryrun: 
			fig.savefig(fname_path, dpi=self.dpi)
		print "Figure created in %s" % fname_path
		plt.close(fig)

	def plot_integrated_correlation_time(self):
		"""
		Plots the integrated correlation through the flowing.
		"""
		# Sets up values to be plotted
		y = self.integrated_autocorrelation_time
		y_std = self.integrated_autocorrelation_time_error
		x = self.a*np.sqrt(8*self.x)

		# Gives title and file name
		title_string = r"$\tau_{int}$ of %s, $N_{cfg}=%2d$" % (self.observable_name, self.N_configurations)
		fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_integrated_ac_time_beta{1:<s}{2:<s}.png".format(self.observable_name_compact, str(self.beta).replace('.','_'), self.fname_addon))

		# Sets up the plot
		fig = plt.figure()
		ax = fig.add_subplot(111)

		# Plot with error
		ax.plot(x,y,color="0")
		ax.fill_between(x, y - y_std, y + y_std, alpha=0.5, edgecolor='', facecolor='#6699ff')
		ax.set_xlabel(r"$\sqrt{8t_{flow}}$")
		ax.set_ylabel(r"$\tau_{int}$")
		ax.set_title(title_string)
		ax.grid(True)
		if not self.dryrun: 
			fig.savefig(fname_path, dpi=self.dpi)
		print "Figure created in %s" % fname_path
		plt.close(fig)

	def plot_histogram(self, flow_time_index, x_label=None, NBins=30, x_limits="equal"):
		"""
		Function for creating histograms of the original, bootstrapped and jackknifed datasets together.
		Args:
			flow_time_index		(int): flow time to plot.
			x_label				(str): x-axis label for plot.
			[optional] NBins	(int): number of histogram bins.
			[optional] x_limits	(str): type of x-axis limits. Default: 'auto'. Choices: 'equal','auto','analysis'

		Raises:
			AssertionError: if the lengths of the different analyzed data sets 
				differ.
			AssertionError: if the flow
		"""
		# Setting proper flow-time 
		if flow_time_index < 0:
			flow_time_index = len(self.unanalyzed_y_data) - abs(flow_time_index)
			assert len(self.unanalyzed_y_data) == len(self.bs_y_data) == len(self.jk_y_data), "Flow lengths of data sets is not equal!"

		# X-label set as the default y-label
		if isinstance(x_label, types.NoneType):
			x_label = self.y_label

		# Ensures flow time is within bounds.
		assert flow_time_index < len(self.unanalyzed_y_data), "Flow time %d is out of bounds." % flow_time_index

		# Sets up title and file name strings
		title_string = r"Spread of %s, $t_{flow}=%.2f$" % (self.observable_name, flow_time_index*self.flow_epsilon)
		fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_histogram_flowt{1:>04d}_beta{2:<s}{3:<s}.png".format(self.observable_name_compact, abs(flow_time_index), str(self.beta).replace('.','_'), self.fname_addon))

		# Sets up plot
		fig = plt.figure()

		# Adds unanalyzed data
		ax1 = fig.add_subplot(311)
		# weights1 = np.ones_like(self.unanalyzed_y_data[flow_time_index])/float(len(self.unanalyzed_y_data[flow_time_index]))
		weights1 = None
		x1, y1, _ = ax1.hist(self.unanalyzed_y_data[flow_time_index], bins=NBins, label="Unanalyzed", weights=weights1)
		ax1.legend()
		ax1.grid("on")
		ax1.set_title(title_string)

		# Adds bootstrapped data
		ax2 = fig.add_subplot(312)
		# weights2 = np.ones_like(self.bs_y_data[flow_time_index])/float(len(self.bs_y_data[flow_time_index]))
		weights2 = None
		x2, y2, _ = ax2.hist(self.bs_y_data[flow_time_index], bins=NBins, label="Bootstrap", weights=weights2)
		ax2.grid("on")
		ax2.legend()
		ax2.set_ylabel("Hits")

		# Adds jackknifed histogram
		ax3 = fig.add_subplot(313)
		# weights3 = np.ones_like(self.jk_y_data[flow_time_index])/float(len(self.jk_y_data[flow_time_index]))
		weights3 = None
		x3, y3, _ = ax3.hist(self.jk_y_data[flow_time_index], bins=NBins, label="Jackknife", weights=weights3)
		ax3.legend()
		ax3.grid("on")
		ax3.set_xlabel(r"%s" % x_label)

		if x_limits == "auto":
			# Lets matplotlib decide on axes
			xlim_positive = None
			xlim_negative = None
		elif x_limits == "equal":
			# Sets the x-axes to be equal
			xlim_positive = np.max([y1, y2, y3])
			xlim_negative = np.min([y1, y2, y3])

			# Sets the axes limits
			ax1.set_xlim(xlim_negative, xlim_positive)
		elif x_limits == "analysis":
			# Sets only the analysises axes equal
			xlim_positive = np.max([y2, y3])
			xlim_negative = np.min([y2, y3])
		else:
			raise KeyError("%s not recognized.\nOptions: 'equal','auto','analysis'." % x_limits)

		# Sets the axes limits
		ax2.set_xlim(xlim_negative, xlim_positive)
		ax3.set_xlim(xlim_negative, xlim_positive)

		# Saves figure
		if not self.dryrun:
			plt.savefig(fname_path, dpi=self.dpi)
		print "Figure created in %s" % fname_path

		# Closes figure for garbage collection
		plt.close(fig)

	def plot_mc_history(self, flow_time_index, correction_function=lambda x: x):
		"""
		Plots the Monte Carlo history at a given flow time .

		Args:
			flow_time_index: index of the flow-time we are to plot for. -1 gives the
				last flow time element.
		
		Raises:
			AssertionError: if the flow_time_index is out of bounds.
		"""

		# Converts flow_time_index if it is minus 1
		if flow_time_index == -1:
			flow_time_index = self.NFlows - 1

		# Ensures flow time is within bounds.
		assert flow_time_index < len(self.unanalyzed_y_data), "Flow time %d is out of bounds." % flow_time_index

		# Sets up title and file name strings
		title_string = r"Monte Carlo history at $t_{flow} = %.2f$" % (flow_time_index*self.flow_epsilon)
		fname_path = os.path.join(self.observable_output_folder_path, "{0:<s}_mchistory_flowt{1:>04d}_beta{2:<s}{3:<s}.png".format(self.observable_name_compact, flow_time_index, str(self.beta).replace('.','_'), self.fname_addon))

		# Sets up plot
		fig = plt.figure()
		ax = fig.add_subplot(111)
		# print self.y_limits
		# ax.set_ylim(self.y_limits)
		ax.plot(correction_function(self.unanalyzed_y_data[flow_time_index]), color="0", label=self.observable_name)
		ax.set_xlabel(r"Monte Carlo time")
		ax.set_ylabel(r"")
		ax.set_title(title_string)
		ax.grid(True)
		ax.legend()
		if not self.dryrun: 
			fig.savefig(fname_path, dpi=self.dpi)
		print "Figure created in %s" % fname_path
		plt.close(fig)

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += "\n" + "="*100
		return return_string


class _AnalyseTopSusBase(FlowAnalyser):
	"""Topological susceptibility analysis base class."""
	observable_name = "Topological Susceptibility"
	observable_name_compact = "topsus"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"

	def __init__(self, *args, **kwargs):
		super(_AnalyseTopSusBase, self).__init__(*args, **kwargs)
		self.__set_size()

	def chi(self, Q_squared):
		"""Topological susceptibility function."""
		return self.const*Q_squared**(0.25)

	def chi_std(self, Q_squared, Q_squared_std):
		"""Topological susceptibility with error propagation."""
		return 0.25*self.const*Q_squared_std / Q_squared**(0.75)

	def __set_size(self):
		"""Function that sets the lattice size deepending on the beta value."""
		# Retrieves lattice spacing
		self.function_derivative = ptools._chi_derivative

		# Sets up constants used in the chi function for topological susceptibility
		self.V = self.lattice_sizes[self.beta]
		self.const = self.hbarc/self.a/self.V**(1./4)
		self.function_derivative_parameters = {"const": self.const}

	def jackknife(self, F=None, F_error=None, store_raw_jk_values=True):
		"""Overriding the jackknife class by adding the chi-function"""
		super(_AnalyseTopSusBase, self).jackknife(F=self.chi,
			F_error=self.chi_std, store_raw_jk_values=store_raw_jk_values)

	def boot(self, N_bs, F=None, F_error=None, store_raw_bs_values=True):
		"""Overriding the bootstrap class by adding the chi-function"""
		super(_AnalyseTopSusBase, self).boot(N_bs, F=self.chi,
			F_error=self.chi_std, store_raw_bs_values=store_raw_bs_values)

class AnalysePlaquette(FlowAnalyser):
	"""Plaquette analysis class."""
	observable_name = "Plaquette"
	observable_name_compact = "plaq"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$P_{\mu\nu}$"

class AnalyseEnergy(FlowAnalyser):
	"""Energy/action density analysis class."""
	observable_name = "Energy"
	observable_name_compact = "energy"
	x_label = r"$t/r_0^2$" # Dimensionsless, Implied multiplication by a^2
	y_label = r"$t^2\langle E \rangle$" # Energy is dimension 4, while t^2 is dimension invsere 4, or length/time which is inverse energy, see Peskin and Schroeder

	def __init__(self, *args, **kwargs):
		super(AnalyseEnergy, self).__init__(*args, **kwargs)
		self.y *= -1.0

		self.x_vals = self.x / self.r0**2 * self.a**2

		self.plot_hline_at = 0.3

	def correction_function(self, y):
		# *  self.flow_epsilon * self.flow_epsilon
		return y * self.x * self.x # factor 0.5 left out, see paper by 

	def plot_original(self):
		super(AnalyseEnergy, self).plot_boot(x=self.x_vals, correction_function=self.correction_function, _plot_bs=False)

	def plot_boot(self):
		super(AnalyseEnergy, self).plot_boot(x=self.x_vals, correction_function=self.correction_function)

	def plot_jackknife(self):
		super(AnalyseEnergy, self).plot_jackknife(x=self.x_vals, correction_function=self.correction_function)

	def plot_histogram(self, *args, **kwargs):
		kwargs["x_label"] = r"$\langle E \rangle$[GeV]"
		super(AnalyseEnergy, self).plot_histogram(*args, **kwargs)

class AnalyseTopologicalSusceptibility(_AnalyseTopSusBase):
	"""Topological susceptibility analysis class."""
	observable_name = "Topological Susceptibility"
	observable_name_compact = "topsus"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"

	def __init__(self, *args, **kwargs):
		super(AnalyseTopologicalSusceptibility, self).__init__(*args, **kwargs)
		self.y **= 2

class AnalyseTopologicalCharge(FlowAnalyser):
	"""Topological charge analysis class."""
	observable_name = "Topological charge"
	observable_name_compact = "topc"
	x_label = r"$\sqrt{8t_{flow}}[fm]$" # Implied multiplication by a
	# y_label = r"$Q = - \sum_x \frac{1}{64 \cdot 32\pi^2}\epsilon_{\mu\nu\rho\sigma}Tr\{G^{clov}_{\mu\nu}G^{clov}_{\rho\sigma}\}$"
	y_label = r"$Q$" # Dimensionsless, a's from discretization cancel out the 1/a^4 in the integration

class AnalyseQQuartic(_AnalyseTopSusBase):
	"""Class for topological susceptibility with quartic topological charge."""
	observable_name = r"$\chi(\langle Q^4 \rangle)^{1/8}$"
	observable_name_compact = "topq4"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	# y_label = r"$\chi(\langle Q^4 \rangle)^{1/8} = \frac{\hbar}{aV^{1/4}} \langle Q^4 \rangle^{1/8} [GeV]$" # 1/8 correct?
	y_label = r"$\chi(\langle Q^4 \rangle)^{1/8} [GeV]$" # 1/8 correct?

	def __init__(self, *args, **kwargs):
		super(AnalyseQQuartic, self).__init__(*args, **kwargs)
		self.y **= 4

	def chi(self, Q4):
		"""Topological susceptibility function."""
		return self.const * Q4**(0.125)

	def chi_std(self, Q4, Q4_std):
		"""Topological susceptibility with error propagation."""
		return 0.125*self.const * Q4_std / Q4**(0.875)

class AnalyseQtQZero(_AnalyseTopSusBase):
	"""Topological susceptibility QtQ0 analysis class."""
	observable_name = r"$\chi(\langle Q_t Q_{t_0} \rangle)^{1/4}$"
	observable_name_compact = "qtqzero"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi(\langle Q_{t} Q_{t_0} \rangle)^{1/4} [GeV]$" # $\chi_t^{1/4}[GeV]$

	def __init__(self, *args, **kwargs):
		super(AnalyseQtQZero, self).__init__(*args, **kwargs)
		self.observable_output_folder_path_old = self.observable_output_folder_path
		self.post_analysis_folder_old = self.post_analysis_folder
		self.y_original = copy.deepcopy(self.y)

	def setQ0(self, q_flow_time_zero_percent, y_label=None):
		"""
		Sets the flow time we are to analyse for
		q_flow_time_zero_percent: float between 0.0 and 1.0, in which we choose what percentage point of the data we set as q0.
		E.g. if it is 0.9, it will be the Q that is closest to 90% of the whole flowed time
		"""

		# Finds the q flow time zero value
		self.q_flow_time_zero = q_flow_time_zero_percent * (self.a * np.sqrt(8*self.x))[-1]
		self.plot_vline_at = self.q_flow_time_zero
		
		# Finds the flow time zero index
		self.flow_time_zero_index = np.argmin(np.abs(self.a * np.sqrt(8*self.x) - self.q_flow_time_zero))

		# Sets file name
		self.observable_name = r"$\chi(\langle Q_t Q_{t_0} \rangle)^{1/4}$ at $t=%.2f$" % (self.q_flow_time_zero)

		# Manual method for multiplying the matrices
		y_q0 = copy.deepcopy(self.y_original[:,self.flow_time_zero_index])
		self.y = copy.deepcopy(self.y_original)

		# Multiplying QtQ0
		for iFlow in xrange(self.y.shape[1]):
			self.y[:,iFlow] *= y_q0

		# self.y = np.abs(self.y)

		if y_label != None:
			self.y_label = y_label

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "%04d" % self.flow_time_zero_index)
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

		# Checks that {post_analysis_folder}/{observable_name}/{time interval} exist
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_old, "%04d" % self.flow_time_zero_index)
		check_folder(self.post_analysis_folder, self.dryrun, self.verbose)

		# Resets some of the ac, jk and bs variable
		self.bootstrap_performed = False
		self.jackknife_performed = False
		self.autocorrelation_performed = False

	def jackknife(self, F=None, F_error=None, store_raw_jk_values=True):
		"""Overriding the jackknife class by adding the Correaltor function"""
		super(AnalyseQtQZero, self).jackknife(F=self.chi,
			F_error=self.chi_std, store_raw_jk_values=store_raw_jk_values)

	def boot(self, N_bs, F=None, F_error=None, store_raw_bs_values=True):
		"""Overriding the bootstrap class by adding the Correaltor function"""
		super(AnalyseQtQZero, self).boot(N_bs, F=self.chi,
			F_error=self.chi_std, store_raw_bs_values=store_raw_bs_values)

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += info_string("Flow time t0", "%.2f" % self.q_flow_time_zero)
		return_string += "\n" + "="*100
		return return_string

class AnalyseQtQ0Euclidean(FlowAnalyser):
	"""Topological susceptibility QtQ0 in euclidean time analysis class."""
	observable_name = r"$\chi(\langle Q_{t_e} Q_{t_{e,0}} \rangle)^{1/4}$"
	observable_name_compact = "qtq0e"
	x_label = r"$t_e[fm]$"
	y_label = r"$\chi(\langle Q_{t_e} Q_{t_{e,0}} \rangle)^{1/4} [GeV]$" # $\chi_t^{1/4}[GeV]$

	def __init__(self, *args, **kwargs):
		super(AnalyseQtQ0Euclidean, self).__init__(*args, **kwargs)
		self.y_original = copy.deepcopy(self.y)
		self.NT = self.y_original.shape[-1]
		
		# Stores old variables for resetting at each new flow time
		self.N_configurations_old = self.N_configurations
		self.observable_output_folder_path_old = self.observable_output_folder_path
		self.post_analysis_folder_old = self.post_analysis_folder
		self.NFlows_old = self.NFlows

	def set_flow_time(self, flow_time_index, euclidean_percent):
		"""Function for setting the flow time we will plot in euclidean time."""

		# Finds the q flow time zero value
		self.flow_time_index = flow_time_index
		self.flow_time = self.x[flow_time_index]
		self.euclidean_time = int((self.NT - 1) * euclidean_percent)

		# Restores y from original data
		self.y = copy.deepcopy(self.y_original)
		self.y = self.y[:,flow_time_index,:]

		# Sets the number of flows as the number of euclidean time slices,
		# as that is now what we are plotting in.
		assert self.y.shape[1] == self.NT, "the first row does not NT."
		self.NFlows = self.NT

		self.V = self.lattice_sizes[self.beta] / float(self.NT)
		self.const = self.hbarc/self.a/self.V**(1./4)
		self.function_derivative_parameters = {"const": self.const}

		# Sets file name
		self.observable_name = r"$\chi(\langle Q_{t_e} Q_{t_{e,0}} \rangle)^{1/4}$ at $t_e=%.2f$, $t_{flow}=%.2f$" % (self.euclidean_time, self.flow_time)
		# self.observable_name = r"$\chi(\langle Q_{t_e} Q_{t_{e,0}} \rangle)^{1/4}$ at $t_{flow}=%.2f$" % (self.flow_time)

		y_e0 = copy.deepcopy(self.y_original[:, self.flow_time_index, self.euclidean_time])
		# print y_e0.shape
		# print self.y.shape

		# Multiplying QtQ0
		for iFlow in xrange(self.NFlows):
			self.y[:,iFlow] *= y_e0

		# print "Taking abs @ 1035"
		# self.y = np.abs(self.y)

		# Sets up variables deependent on the number of configurations again
		self.unanalyzed_y_data = np.zeros((self.NFlows, self.N_configurations))
		self.autocorrelations = np.zeros((self.NFlows, self.N_configurations/2))
		self.autocorrelations_errors = np.zeros((self.NFlows, self.N_configurations/2))

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "%04d" % self.flow_time_index)
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

		# Checks that {post_analysis_folder}/{observable_name}/{time interval} exist
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_old, "%04d" % self.flow_time_index)
		check_folder(self.post_analysis_folder, self.dryrun, self.verbose)

		# Resets some of the ac, jk and bs variable
		self.bootstrap_performed = False
		self.jackknife_performed = False
		self.autocorrelation_performed = False

	def C(self, qtq0):
		"""Correlator for qtq0."""
		return self.const*qtq0

	def C_std(self, qtq0, qtq0_std):
		"""Correlator for qtq0 with error propagation."""
		return self.const*qtq0_std

	def jackknife(self, F=None, F_error=None, store_raw_jk_values=True):
		"""Overriding the jackknife class by adding the Correaltor function"""
		super(AnalyseQtQ0Euclidean, self).jackknife(F=self.C,
			F_error=self.C_std, store_raw_jk_values=store_raw_jk_values)

	def boot(self, N_bs, F=None, F_error=None, store_raw_bs_values=True):
		"""Overriding the bootstrap class by adding the Correaltor function"""
		super(AnalyseQtQ0Euclidean, self).boot(N_bs, F=self.C,
			F_error=self.C_std, store_raw_bs_values=store_raw_bs_values)

	# def autocorrelation(self, store_raw_ac_error_correction=True, method="wolff"):
	# 	"""Overriding the ac class."""
	# 	super(AnalyseQtQ0Euclidean, self).autocorrelation(store_raw_ac_error_correction=True, method="luscher")

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += info_string("Flow time:", "%f" % self.flow_time)
		return_string += "\n" + "="*100
		return return_string


class AnalyseTopChargeInEuclTime(FlowAnalyser):
	"""Analysis of the topological charge in Euclidean Time."""
	observable_name = "Topological Charge in Euclidean Time"
	observable_name_compact = "topct"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q_{t_{euclidean}}$"

	def __init__(self, *args, **kwargs):
		super(AnalyseTopChargeInEuclTime, self).__init__(*args, **kwargs)
		self.y_original = copy.deepcopy(self.y)
		self.observable_output_folder_path_old = self.observable_output_folder_path
		self.post_analysis_folder_old = self.post_analysis_folder
		self.NT = self.y_original.shape[-1]

	def setEQ0(self, t_euclidean_index):
		"""
		Sets the Euclidean time we are to analyse for. Q_{t_E=}
		q_flow_time_zero_percent: float between 0.0 and 1.0, in which we choose what percentage point of the data we set as q0.
		E.g. if it is 0.9, it will be the Q that is closest to 90% of the whole flowed time.

		Args:
			t_euclidean_index: integer of what time point we will look at
		"""
		# Finds the euclidean time zero index
		self.t_euclidean_index = t_euclidean_index

		# Sets file name
		self.observable_name = r"$Q_{t_{euclidean}}$ at $i_{euclidean}=%d$" % self.t_euclidean_index

		# Manual method for multiplying the matrices
		self.y = copy.deepcopy(self.y_original[:,:,self.t_euclidean_index])

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "%04d" % self.t_euclidean_index)
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

		# Checks that {post_analysis_folder}/{observable_name}/{time interval} exist
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_old, "%04d" % self.t_euclidean_index)
		check_folder(self.post_analysis_folder, self.dryrun, self.verbose)

		# Resets some of the ac, jk and bs variable
		self.bootstrap_performed = False
		self.jackknife_performed = False
		self.autocorrelation_performed = False

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += info_string("Euclidean time", "%d" % self.t_euclidean_index)
		return_string += "\n" + "="*100
		return return_string


class AnalyseTopSusInEuclTime(_AnalyseTopSusBase):
	"""Analysis of the topological susceptibility in Euclidean Time."""
	observable_name = "Topological Susceptibility in Euclidean Time"
	observable_name_compact = "topsust"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi(\langle Q_t Q_{t_{euclidean}} \rangle)^{1/4} [GeV]$"

	def __init__(self, *args, **kwargs):
		super(AnalyseTopSusInEuclTime, self).__init__(*args, **kwargs)
		self.y_original = copy.deepcopy(self.y)
		self.observable_output_folder_path_old = self.observable_output_folder_path
		self.NT = self.y_original.shape[-1]

	def setEQ0(self, t_euclidean_index):
		"""
		Sets the Euclidean time we are to analyse for. Q_{t_E=}
		q_flow_time_zero_percent: float between 0.0 and 1.0, in which we choose what percentage point of the data we set as q0.
		E.g. if it is 0.9, it will be the Q that is closest to 90% of the whole flowed time.

		Args:
			t_euclidean_index: integer of what time point we will look at
		"""

		# Finds the euclidean time zero index
		self.t_euclidean_index = t_euclidean_index

		self.V = self.lattice_sizes[self.beta] / float(self.NT)
		self.const = self.hbarc/self.a/self.V**(1./4)
		self.function_derivative_parameters = {"const": self.const}

		# Sets file name
		self.observable_name = r"$\chi(\langle Q_t Q_{t_{euclidean}} \rangle)^{1/4}$ at $i_{euclidean}=%d$" % self.t_euclidean_index

		# Manual method for multiplying the matrices
		y_qe0 = copy.deepcopy(self.y_original[:,:,self.t_euclidean_index])
		self.y = copy.deepcopy(self.y_original)

		# Sums the euclidean time
		self.y = np.sum(self.y, axis=2)

		self.y *= y_qe0

		# self.y = np.abs(self.y) # If not absolute value, will get error!

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "%04d" % self.t_euclidean_index)
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

		# Checks that {post_analysis_folder}/{observable_name}/{time interval} exist
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_old, "%04d" % self.t_euclidean_index)
		check_folder(self.post_analysis_folder, self.dryrun, self.verbose)

		# Resets some of the ac, jk and bs variable
		self.bootstrap_performed = False
		self.jackknife_performed = False
		self.autocorrelation_performed = False

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += info_string("Euclidean time", "%d" % self.t_euclidean_index)
		return_string += "\n" + "="*100
		return return_string

	def plot_histogram(self, *args, **kwargs):
		print "Not plotting histogram for %s at euclidean time %d due to weird error." % (self.observable_name_compact, self.t_euclidean_index)

class _EuclideanSplitAnalysis(FlowAnalyser):
	"""
	Base class for splitting in either Euclidean and creating separate data
	sets from that.
	"""
	def __init__(self, *args, **kwargs):
		super(_EuclideanSplitAnalysis, self).__init__(*args, **kwargs)
		self.y_original = copy.deepcopy(self.y)

	def set_t_interval(self, t_interval):
		self.y = copy.deepcopy(self.y_original)
		self.y = self.y[:,:,t_interval[0]:t_interval[1]]
		self.t_interval = t_interval

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "int%03d-%03d" % t_interval)
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

		# Checks that {post_analysis_folder}/{observable_name}/{time interval} exist
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_old, "%03d-%03d" % self.t_interval)
		check_folder(self.post_analysis_folder, self.dryrun, self.verbose)

		# Resets some of the ac, jk and bs variable
		self.bootstrap_performed = False
		self.jackknife_performed = False
		self.autocorrelation_performed = False

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += info_string("Time interval", "[%d,%d)" % self.t_interval)
		return_string += "\n" + "="*100
		return return_string


class AnalyseTopChargeSplitEuclideanTime(_EuclideanSplitAnalysis):
	"""
	Analysis where one can split the topological in Euclidean time to obtain 
	an estimate of the topological charge.
	"""
	observable_name = "Topological Charge in Euclidean Time"
	observable_name_compact = "topcte"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q$"

	def __init__(self, *args, **kwargs):
		super(AnalyseTopChargeSplitEuclideanTime, self).__init__(*args, **kwargs)
		self.NT = self.y_original.shape[-1]

		self.post_analysis_folder_old = self.post_analysis_folder
		self.observable_output_folder_path_old = self.observable_output_folder_path

	def set_t_interval(self, *args):
		"""Runs first the inherited time setter function, then its own."""
		super(AnalyseTopChargeSplitEuclideanTime, self).set_t_interval(*args)
		self.y = np.sum(self.y, axis=2)
		self.observable_name = r"Q in Euclidean time $t=[%d,%d)$" % self.t_interval

class AnalysisTopSusSplitEuclideanTime(_EuclideanSplitAnalysis, _AnalyseTopSusBase):
	"""
	Analysis where one can split the topological susceptibility in Euclidean
	time to obtain an estimate of the topological charge.
	"""
	observable_name = "Topological Susceptibility in MC Time"
	observable_name_compact = "topsuste"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi^{1/4} [GeV]$"

	def __init__(self, *args, **kwargs):
		super(AnalysisTopSusSplitEuclideanTime, self).__init__(*args, **kwargs)
		self.NT = self.y_original.shape[-1]
		self.observable_output_folder_path_old = self.observable_output_folder_path

	def set_t_interval(self, *args):
		"""Runs first the inherited time setter function, then its own."""
		super(AnalysisTopSusSplitEuclideanTime, self).set_t_interval(*args)
		self.observable_name = r"$\chi(\langle Q^2 \rangle)^{1/4}$ in Euclidean time $[%d,%d)$" % self.t_interval
		self.NT_interval_size = self.t_interval[-1] - self.t_interval[0]
		self.V = self.lattice_sizes[self.beta] * self.NT_interval_size / float(self.NT)
		self.const = self.hbarc/self.a/self.V**(1./4)
		self.function_derivative_parameters = {"const": self.const}

		self.y = np.sum(self.y, axis=2)
		self.y **= 2

class _MCSplitAnalysis(FlowAnalyser):
	"""Base class for analysing different parts of the Monte Carlo history."""
	def __init__(self, *args, **kwargs):
		super(_MCSplitAnalysis, self).__init__(*args, **kwargs)
		self.y_original = copy.deepcopy(self.y)
		self.N_configurations_old = self.N_configurations

	def set_MC_interval(self, MC_interval):
		self.y = copy.deepcopy(self.y_original)
		self.y = self.y[MC_interval[0]:MC_interval[1],:]

		self.N_configurations = self.y.shape[0]

		# Sets up variables deependent on the number of configurations again
		self.unanalyzed_y_data = np.zeros((self.NFlows, self.N_configurations))
		self.autocorrelations = np.zeros((self.NFlows, self.N_configurations/2))
		self.autocorrelations_errors = np.zeros((self.NFlows, self.N_configurations/2))

		self.MC_interval = MC_interval

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "MCint%03d-%03d" % MC_interval)
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

		# Checks that {post_analysis_folder}/{observable_name}/{time interval} exist
		self.post_analysis_folder = os.path.join(self.post_analysis_folder_old, "%03d-%03d" % self.MC_interval)
		check_folder(self.post_analysis_folder, self.dryrun, self.verbose)

		# Resets some of the ac, jk and bs variable
		self.bootstrap_performed = False
		self.jackknife_performed = False
		self.autocorrelation_performed = False

	def __str__(self):
		info_string = lambda s1, s2: "\n{0:<20s}: {1:<20s}".format(s1, s2)
		return_string = ""
		return_string += "\n" + "="*100
		return_string += info_string("Data batch folder", self.batch_data_folder)
		return_string += info_string("Batch name", self.batch_name)
		return_string += info_string("Observable", self.observable_name_compact)
		return_string += info_string("Beta", "%.2f" % self.beta)
		return_string += info_string("Time interval", "[%d,%d)" % self.MC_interval)
		return_string += "\n" + "="*100
		return return_string


class AnalysisTopChargeSplitMCTime(_MCSplitAnalysis):
	"""
	Analysis where one can split the topological charge in Monte Carlo time to
	obtain an estimate of the topological charge.
	"""
	observable_name = "Topological Charge in MC Time"
	observable_name_compact = "topcMC"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q$"

	def __init__(self, *args, **kwargs):
		super(AnalysisTopChargeSplitMCTime, self).__init__(*args, **kwargs)
		self.NT = self.y_original.shape[-1]
		self.observable_output_folder_path_old = self.observable_output_folder_path

	def set_MC_interval(self, *args):
		"""Runs first the inherited time setter function, then its own."""
		super(AnalysisTopChargeSplitMCTime, self).set_MC_interval(*args)
		self.observable_name = r"Q in MC interval $[%d,%d)$" % self.MC_interval

class AnalysisTopSusSplitMCTime(_MCSplitAnalysis, _AnalyseTopSusBase):
	"""
	Analysis where one can split the topological susceptibility in Monte Carlo
	time to obtain an estimate of the topological charge.
	"""
	observable_name = "Topological Susceptibility in MC Time"
	observable_name_compact = "topsusMC"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi^{1/4} [GeV]$"

	def __init__(self, *args, **kwargs):
		super(AnalysisTopSusSplitMCTime, self).__init__(*args, **kwargs)
		self.NT = self.y_original.shape[-1]

		# Squares the values
		self.y_original **= 2

		self.observable_output_folder_path_old = self.observable_output_folder_path

	def set_MC_interval(self, *args):
		"""Runs first the inherited time setter function, then its own."""
		super(AnalysisTopSusSplitMCTime, self).set_MC_interval(*args)
		self.observable_name = r"$\chi(\langle Q^2 \rangle)^{1/4}$ in MC interval $[%d,%d)$" % self.MC_interval


def main():
	print "Module FlowAnalyser intended to be inherited, and not used as a standalone analysis tool."
	exit(1)

if __name__ == '__main__':
	main()