from tools.folderreadingtools import GetDirectoryTree, GetFolderContents, write_data_to_file, check_folder, write_raw_analysis_to_file
import statistics.parallel_tools as ptools
from statistics.jackknife import Jackknife
from statistics.bootstrap import Bootstrap
from statistics.autocorrelation import Autocorrelation
import os, numpy as np, matplotlib.pyplot as plt, sys, pandas as pd, multiprocessing, pickle, time, copy

# from tqdm import trange

class FlowAnalyser(object):
	observable_name = "Missing_Observable_Name"
	observable_name_compact = "missing_obs_name"
	x_label = "Missing x-label"
	y_label = "Missing y-label"
	mark_interval = 5
	error_mark_interval = 5
	autocorrelations_limits = 1
	figures_folder = "../figures"
	fname_addon = ""
	function_derivative = None
	dpi = 350

	def __init__(self, file_tree, batch_name, data=None, dryrun=False, flow=True, parallel = False, numprocs = 4, verbose = False, figures_folder = False, create_perflow_data = False):
		# Sets up global constants
		self.batch_name = batch_name
		self.batch_data_folder = file_tree.data_batch_folder
		self.N_bs = None
		self.dryrun = dryrun
		self.flow = flow
		self.verbose = verbose
		if figures_folder != False:
			self.figures_folder = figures_folder

		# Parallel variables
		self.parallel = parallel
		self.numprocs = numprocs

		# Checks that a figures folder exists
		check_folder(self.figures_folder,self.dryrun,verbose=self.verbose)

		# Check that a data run folder exist, so one data anlysis performed on different data sets do not mix
		self.data_batch_folder_path = os.path.join(self.figures_folder,os.path.split(self.batch_data_folder)[-1])
		check_folder(self.data_batch_folder_path,self.dryrun,verbose=self.verbose)

		# Checks that a batch folder exists
		self.batch_name_folder_path = os.path.join(self.data_batch_folder_path,self.batch_name)
		check_folder(self.batch_name_folder_path,self.dryrun,verbose=self.verbose)

		# Checks that observable output folder exist, and if not will create it
		self.observable_output_folder_path = os.path.join(self.batch_name_folder_path,self.observable_name_compact)
		check_folder(self.observable_output_folder_path,self.dryrun,verbose=self.verbose)

		# Prints observable and batch
		print "="*100
		print "Batch:       %s" % self.batch_name
		print "Observables: %s" % self.observable_name_compact
		print "="*100

		# Enables possibility of providing data
		if data == None:
			self.data = GetFolderContents(file_tree,self.observable_name_compact,flow=self.flow)
		else:
			self.data = copy.deepcopy(data)

		# Creates perflow data if prompted
		if create_perflow_data:
			self.data.create_perflow_data(verbose=self.verbose)

		if self.verbose:
			print "Size of data: %.4g kB" % (sys.getsizeof(self.data.data_y)/1024.0)

		# Sets up variables
		self.y = self.data.data_y
		self.x = self.data.data_x

		# Max plotting window variables
		self.y_limits = [None,None]

		self.N_configurations, self.NFlows = self.y.shape

		# Small error checking in retrieving number of flows
		if (int(self.data.meta_data["NFlows"]) + 1) != self.NFlows:
			raise ValueError("Number of flows %d does not match the number provided by metadata %d." % (self.NFlows,int(self.data.meta_data["NFlows"]) + 1))

		# Non-bootstrapped data
		self.unanalyzed_y = np.zeros(self.NFlows)
		self.unanalyzed_y_std = np.zeros(self.NFlows)
		self.unanalyzed_y_data = np.zeros((self.NFlows,self.N_configurations))

		# Bootstrap data
		self.bootstrap_performed = False
		self.bs_y = np.zeros(self.NFlows)
		self.bs_y_std = np.zeros(self.NFlows)

		# Jackknifed data
		self.jackknife_performed = False
		self.jk_y = np.zeros(self.NFlows)
		self.jk_y_std = np.zeros(self.NFlows)
		self.jk_y_data = np.zeros((self.NFlows,self.N_configurations))

		# Autocorrelation data
		self.autocorrelation_performed = False
		self.autocorrelations = np.zeros((self.NFlows,self.N_configurations/2))
		self.autocorrelations_errors = np.zeros((self.NFlows,self.N_configurations/2))
		self.integrated_autocorrelation_time = np.ones(self.NFlows)
		self.integrated_autocorrelation_time_error = np.zeros(self.NFlows)
		self.autocorrelation_error_correction = np.ones(self.NFlows)

		# Gets the lattice spacing
		self.beta = self.data.meta_data["beta"]
		self.a = self.getLatticeSpacing(self.beta)
		self.r = 0.5 # Sommer Parameters

	def __check_ac(self,fname):
		"""
		If autocorrelation has been performed it will add "_noErrorCorrection" to the filename
		Args:
			filename (str)	: figure file name
		Returns:
			filename (str)	: possibly modified figure file name
		"""
		head, ext = os.path.splitext(fname)
		fname_addon = "_noErrorCorrection"
		if not self.autocorrelation_performed:
			return head + "_noErrorCorrection" + ext
		else:
			return head + ext
			
	@staticmethod
	def getLatticeSpacing(beta):
		"""
		Static method for calculating the lattice spacing a from the beta-value. Based on paper by M. Guagnelli et al(1998)
		Args:
			beta (float)		: beta value
		Returns:
			a (float)			: lattice spacing
		"""
		if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
		r = 0.5
		bval = (beta - 6)
		a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
		return a

	def boot(self, N_bs, F = ptools._default_return, F_error = ptools._default_error_return, store_raw_bs_values = True):
		# Stores number of bootstraps
		self.N_bs = N_bs
		
		# Sets up raw bootstrap values
		self.bs_y_data = np.zeros((self.NFlows,self.N_bs))

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
			results = pool.map(ptools._bootstrap_parallel_core,input_values)

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

			# Non-parallel method for calculating jackknife
			for i in xrange(self.NFlows):
				bs = Bootstrap(self.y[:,i], N_bs, index_lists = index_lists)
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
		self.unanalyzed_y_std = F_error(self.unanalyzed_y,self.unanalyzed_y_std)
		self.unanalyzed_y = F(self.unanalyzed_y)

		# Runs bs data through function F
		self.bs_y_data = F(self.bs_y_data)
		self.unanalyzed_y_data = F(self.unanalyzed_y_data)

		if store_raw_bs_values:
			# Store as binary
			write_raw_analysis_to_file(self.bs_y_data, "bootstrap", self.observable_name_compact, self.batch_data_folder, self.beta, dryrun = self.dryrun)

		# Sets performed flag to true
		self.bootstrap_performed = True

	def jackknife(self, F = ptools._default_return, F_error = ptools._default_error_return, store_raw_jk_values=True):
		if self.parallel:
			# Sets up jobs for parallel processing
			input_values = [self.y[:,i] for i in xrange(self.NFlows)]

			# Prints job size
			# if self.verbose:
			# 	print sys.getsizeof(input_values)/1024.0, "kB"

			# Initializes multiprocessing
			pool = multiprocessing.Pool(processes=self.numprocs)								

			# Runs parallel processes. Can this be done more efficiently?
			results = pool.map(ptools._jackknife_parallel_core,input_values)

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
				jk = Jackknife(self.y[:,i], jk_statistics = jk_statistics, non_jk_statistics = non_jk_statistics)
				self.jk_y[i] = jk.jk_avg
				self.jk_y_std[i] = jk.jk_std
				self.jk_y_data[i] = jk.jk_data

		# Runs data through the F and F_error
		self.jk_y_std = F_error(self.jk_y,self.jk_y_std)
		self.jk_y = F(self.jk_y)
		self.jk_y_data = F(self.jk_y_data)

		if store_raw_jk_values:
			# Store as binary
			write_raw_analysis_to_file(self.jk_y_data, "jackknife", self.observable_name_compact, self.batch_data_folder, self.beta, dryrun = self.dryrun)

		# Sets performed flag to true
		self.jackknife_performed = True

	def autocorrelation(self, store_raw_ac_error_correction = True):
		"""
		Function for running the autocorrelation routine.
		"""
		# Gets autocorrelation
		if self.parallel:
			# Sets up parallel job
			pool = multiprocessing.Pool(processes=self.numprocs)
			
			if self.function_derivative != None:
				# Sets up jobs for parallel processing
				input_values = zip(	[self.y[:,i] for i in xrange(self.NFlows)],
									[self.function_derivative for i in xrange(self.NFlows)])

				# Initiates parallel jobs
				results = pool.map(ptools._autocorrelation_propagated_parallel_core,input_values)
			else:
				# Sets up jobs for parallel processing
				input_values = zip([self.y[:,i] for i in xrange(self.NFlows)])
				
				# Initiates parallel jobs
				results = pool.map(ptools._autocorrelation_parallel_core,input_values)

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
				if self.function_derivative != None:
					ac = PropagatedAutocorrelation(self.y[:,i])
				else:
					ac = Autocorrelation(self.y[:,i])
				self.autocorrelations[i] = ac.R
				self.autocorrelations_errors[i] = ac.R_error
				self.integrated_autocorrelation_time[i] = ac.integrated_autocorrelation_time()
				self.integrated_autocorrelation_time_error[i] = ac.integrated_autocorrelation_time_error()
				self.autocorrelation_error_correction[i] = np.sqrt(2*self.integrated_autocorrelation_time[i])

				if (i % 10 == 0):
					# Small progressbar
					sys.stdout.write("\rCalculating autocorrelation: %4.1f%% done" % (100*float(i)/float(self.NFlows)))
					sys.stdout.flush()

			# Finalizes
			sys.stdout.write("\rCalculating autocorrelation: 100.0%% done")
			sys.stdout.flush()

		# Stores the ac error correction
		if store_raw_ac_error_correction:
			write_raw_analysis_to_file(self.autocorrelation_error_correction,"autocorrelation",self.observable_name_compact,self.batch_data_folder,self.beta,dryrun=self.dryrun)

		# Sets performed flag to true
		self.autocorrelation_performed = True

	def plot_jackknife(self, x = None, correction_function = lambda x : x):
		# Checks that jacknifing has been performed.
		if not self.jackknife_performed:
			raise ValueError("Jackknifing has not been performed yet.")

		# Sets up the x axis array to be plotted
		if type(x) != np.ndarray:
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x*self.data.meta_data["FlowEpsilon"])
		
		# Copies data over to arrays to be plotted
		y = self.jk_y
		y_std = self.jk_y_std*self.autocorrelation_error_correction

		# Sets up the title and filename strings
		title_string = r"Jacknife of %s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta)
		fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_jackknife_beta{1:<s}{2:<s}.png".format(self.observable_name_compact,str(self.beta).replace('.','_'),self.fname_addon))

		# Plots the jackknifed data
		self.__plot_error_core(x,correction_function(y),correction_function(y_std),title_string,fname_path)

	def plot_boot(self,plot_bs=True, x = None, correction_function = lambda x : x):
		# Checks that the bootstrap has been performed.
		if not self.bootstrap_performed and plot_bs:
			raise ValueError("Bootstrap has not been performed yet.")

		# Retrieves relevant data and sets up the arrays to be plotted
		if type(x) != np.ndarray:
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x*self.data.meta_data["FlowEpsilon"])
		
		# Determines if we are to plot bootstrap or original and retrieves data
		if plot_bs:
			y = self.bs_y
			y_std = self.bs_y_std*self.autocorrelation_error_correction
		else:
			y = self.unanalyzed_y
			y_std = self.unanalyzed_y_std*self.autocorrelation_error_correction

		# Sets up the title and filename strings
		if plot_bs:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta,self.N_bs)
			fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_bootstrap_Nbs{2:<d}_beta{1:<s}{3:<s}.png".format(self.observable_name_compact,str(self.beta).replace('.','_'),self.N_bs,self.fname_addon))
		else:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name, self.data.meta_data["NFlows"],self.beta)
			fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_original_beta{1:<s}{2:<s}.png".format(self.observable_name_compact,str(self.beta).replace('.','_'),self.fname_addon))

		# Plots either bootstrapped or regular stuff
		self.__plot_error_core(x,correction_function(y),correction_function(y_std),title_string,fname_path)

	def plot_original(self, x = None, correction_function = lambda x : x):
		"""
		Plots the default analysis, mean and std of the observable.
		"""
		self.plot_boot(plot_bs=False, x = x, correction_function = correction_function)

	def __plot_error_core(self,x,y,y_std,title_string,fname):
		# Plots the jackknifed data
		fig = plt.figure()
		ax = fig.add_subplot(111)

		# Plots the error bar
		ax.errorbar(x,y,yerr=y_std,fmt=".",color="0",ecolor="r",label=self.observable_name,markevery=self.mark_interval,errorevery=self.error_mark_interval)

		# if self.y_limits[0] == None and self.y_limits[1] == None:
		# 	self.y_limits = [-np.max(y),np.max(y)]

		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)
		ax.set_ylim(self.y_limits)
		ax.grid(True)
		ax.set_title(title_string)
		if not self.dryrun: 
			fig.savefig(self.__check_ac(fname),dpi=self.dpi)
		if self.verbose:
			print "Figure created in %s" % fname
		plt.close(fig)

	def plot_autocorrelation(self,flow_time,plot_abs_value=False):
		# Checks that autocorrelations has been performed.
		if not self.autocorrelation_performed:
			raise ValueError("Autocorrelation has not been performed yet.")

		# sets up the autocorrelation
		N_autocorr = self.N_configurations / 2

		# Converts flow_time if it is minus 1
		if flow_time == -1:
			flow_time = self.NFlows - 1

		# Ensures flow time is within bounds.
		assert flow_time < self.NFlows, "Flow time %d is out of bounds." % flow_time

		# Finds the maximum value at each MC time and sets up the y array
		x = range(N_autocorr)
		y = np.zeros(N_autocorr)
		if plot_abs_value:
			y = np.abs(self.autocorrelations[flow_time,:])
		else:
			y = self.autocorrelations[flow_time,:]
		
		y_std = self.autocorrelations_errors[flow_time,:]

		# Sets up the title and filename strings
		title_string = r"Autocorrelation of %s at flow time $t=%.2f$, $\beta=%.2f$, $N_{cfg}=%2d$" % (self.observable_name,flow_time*self.data.meta_data["FlowEpsilon"],self.beta,self.N_configurations)
		fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_autocorrelation_flowt{1:<d}_beta{2:<s}{3:<s}.png".format(self.observable_name_compact,flow_time, str(self.beta).replace('.','_'),self.fname_addon))

		# Plots the autocorrelations
		fig = plt.figure()
		ax = fig.add_subplot(111)
		# ax.plot(x,y,color="0",label=self.observable_name)
		ax.errorbar(x,y,yerr=y_std,color="0",ecolor="r")#,label=self.observable_name)
		ax.set_ylim(-self.autocorrelations_limits,self.autocorrelations_limits)
		ax.set_xlim(0,N_autocorr)
		ax.set_xlabel(r"Lag $h$")
		ax.set_ylabel(r"$R = \frac{C_h}{C_0}$")
		ax.set_title(title_string)
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.2))
		ax.grid(True)
		# ax.legend()
		if not self.dryrun: 
			fig.savefig(fname_path,dpi=self.dpi)
		print "Figure created in %s" % fname_path
		plt.close(fig)

	def plot_integrated_correlation_time(self):
		"""
		Plots the integrated correlation through the flowing.
		"""
		# Sets up values to be plotted
		y = self.integrated_autocorrelation_time
		y_std = self.integrated_autocorrelation_time_error
		x = self.x*self.data.meta_data["FlowEpsilon"]

		# Gives title and file name
		title_string = r"Integrated autocorrelation time of %s for $\beta=%.2f$, $N_{cfg}=%2d$" % (self.observable_name,self.beta,self.N_configurations)
		fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_integrated_ac_time_beta{1:<s}{2:<s}.png".format(self.observable_name_compact, str(self.beta).replace('.','_'),self.fname_addon))

		# Sets up the plot
		fig = plt.figure()
		ax = fig.add_subplot(111)

		# Plot with error
		ax.plot(x,y,color="0")
		ax.fill_between(x, y - y_std, y + y_std,alpha=0.5, edgecolor='', facecolor='#6699ff')

		# ax.plot(x,y,color="0",label=self.observable_name)
		# ax.errorbar(x,y,yerr=y_std,color="0",ecolor="r")#,label=self.observable_name)
		ax.set_xlabel(r"$t_{flow}$")
		ax.set_ylabel(r"$\tau_{int}$")
		ax.set_title(title_string)
		ax.grid(True)
		if not self.dryrun: 
			fig.savefig(fname_path,dpi=self.dpi)
		print "Figure created in %s" % fname_path
		plt.close(fig)

	def plot_histogram(self, flow_time, x_label, Nbins = 30, x_limits="auto"):
		"""
		Function for creating histograms of the original, bootstrapped and jackknifed datasets together.
		Args:
			flow_time			(int): flow time to plot.
			x_label				(str): x-axis label for plot.
			[optional] Nbins	(int): number of histogram bins.
			[optional] x_limits	(str): type of x-axis limits. Default: 'auto'. Choices: 'equal','auto','analysis'
		"""
		# Setting proper flow-time 
		if flow_time < 0:
			flow_time = len(self.unanalyzed_y_data) - abs(flow_time)
			assert len(self.unanalyzed_y_data) == len(self.bs_y_data) == len(self.jk_y_data), "Flow lengths of data sets is not equal!"
		
		# Ensures flow time is within bounds.
		assert flow_time < len(self.unanalyzed_y_data), "Flow time %d is out of bounds." % flow_time
		
		# Sets up title and file name strings
		title_string = r"Spread of %s, $\beta=%.2f$, flow time $t=%.2f$" % (self.observable_name, float(self.data.meta_data["beta"]),flow_time*self.data.meta_data["FlowEpsilon"])
		fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_histogram_flowt{1:<d}_beta{2:<s}{3:<s}.png".format(self.observable_name_compact,abs(flow_time),str(self.beta).replace('.','_'),self.fname_addon))

		# Sets up plot
		fig = plt.figure()

		# Adds unanalyzed data
		ax1 = fig.add_subplot(311)
		# weights1 = np.ones_like(self.unanalyzed_y_data[flow_time])/float(len(self.unanalyzed_y_data[flow_time]))
		weights1 = None
		x1, y1, _ = ax1.hist(self.unanalyzed_y_data[flow_time],bins=Nbins,label="Unanalyzed", weights = weights1)
		ax1.legend()
		ax1.grid("on")
		ax1.set_title(title_string)

		# Adds bootstrapped data
		ax2 = fig.add_subplot(312)
		# weights2 = np.ones_like(self.bs_y_data[flow_time])/float(len(self.bs_y_data[flow_time]))
		weights2 = None
		x2, y2, _ = ax2.hist(self.bs_y_data[flow_time],bins=Nbins,label="Bootstrap", weights = weights2)
		ax2.grid("on")
		ax2.legend()
		ax2.set_ylabel("Hits")

		# Adds jackknifed histogram
		ax3 = fig.add_subplot(313)
		# weights3 = np.ones_like(self.jk_y_data[flow_time])/float(len(self.jk_y_data[flow_time]))
		weights3 = None
		x3, y3, _ = ax3.hist(self.jk_y_data[flow_time],bins=Nbins,label="Jackknife", weights = weights3)
		ax3.legend()
		ax3.grid("on")
		ax3.set_xlabel(r"%s" % x_label)

		if x_limits == "auto":
			# Lets matplotlib decide on axes
			xlim_positive = None
			xlim_negative = None
		elif x_limits == "equal":
			# Sets the x-axes to be equal
			xlim_positive = np.max([np.max(_y) for _y in [np.abs(y1),np.abs(y2),np.abs(y3)]])
			xlim_negative = -xlim_positive
			
			# Sets the axes limits
			ax1.set_xlim(xlim_negative,xlim_positive)
		elif x_limits == "analysis":
			# Sets only the analysises axes equal
			xlim_positive = np.max([np.max(_y) for _y in [np.abs(y2),np.abs(y3)]])
			xlim_negative = -xlim_positive
		else:
			raise KeyError("%s not recognized.\nOptions: 'equal','auto','analysis'." % x_limits)

		# Sets the axes limits
		ax2.set_xlim(xlim_negative,xlim_positive)
		ax3.set_xlim(xlim_negative,xlim_positive)

		# Saves figure
		if not self.dryrun:
			plt.savefig(fname_path)
		print "Figure created in %s" % fname_path

		# Closes figure for garbage collection
		plt.close(fig)

	def plot_mc_history(self,flow_time,correction_function = lambda x : x):
		"""
		Plots the Monte Carlo history at a given flow time 
		"""
		# Converts flow_time if it is minus 1
		if flow_time == -1:
			flow_time = self.NFlows - 1

		# Ensures flow time is within bounds.
		assert flow_time < len(self.unanalyzed_y_data), "Flow time %d is out of bounds." % flow_time

		# Sets up title and file name strings
		title_string = r"Monte Carlo history for %s, $\beta=%.2f$, $t_{flow} = %.2f$" % (self.observable_name,self.beta,flow_time*self.data.meta_data["FlowEpsilon"])
		fname_path = os.path.join(self.observable_output_folder_path,"{0:<s}_mchistory_flowt{1:<d}_beta{2:<s}{3:<s}.png".format(self.observable_name_compact,flow_time,str(self.beta).replace('.','_'),self.fname_addon))

		# Sets up plot
		fig = plt.figure()
		ax = fig.add_subplot(111)
		# print self.y_limits
		# ax.set_ylim(self.y_limits)
		ax.plot(correction_function(self.unanalyzed_y_data[flow_time]),color="0",label=self.observable_name)
		ax.set_xlabel(r"Monte Carlo time")
		ax.set_ylabel(r"")
		ax.set_title(title_string)
		ax.grid(True)
		ax.legend()
		if not self.dryrun: 
			fig.savefig(fname_path,dpi=self.dpi)
		print "Figure created in %s" % fname_path
		plt.close(fig)

class AnalysePlaquette(FlowAnalyser):
	"""
	Plaquette analysis class.
	"""
	observable_name = "Plaquette"
	observable_name_compact = "plaq"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$P_{\mu\nu}$"

class AnalyseTopologicalCharge(FlowAnalyser):
	"""
	Topological charge analysis class.
	"""
	observable_name = "Topological charge"
	observable_name_compact = "topc"
	x_label = r"$\sqrt{8t_{flow}}[fm]$" # Implied multiplication by a
	y_label = r"$Q = \sum_x \frac{1}{32\pi^2}\epsilon_{\mu\nu\rho\sigma}Tr\{G^{clov}_{\mu\nu}G^{clov}_{\rho\sigma}\}$[GeV]"

class AnalyseEnergy(FlowAnalyser):
	"""
	Energy/action density analysis class.
	"""
	observable_name = "Energy"
	observable_name_compact = "energy"
	x_label = r"$t/r_0^2$" # Dimensionsless, Implied multiplication by a^2
	y_label = r"$t^2\langle E \rangle$" # Energy is dimension 4, while t^2 is dimension invsere 4, or length/time which is inverse energy, see Peskin and Schroeder

	def __init__(self,*args,**kwargs):
		super(AnalyseEnergy,self).__init__(*args,**kwargs)
		self.y *= -self.y/64.0

	def correction_function(self, y):
		return y*self.x*self.x*self.data.meta_data["FlowEpsilon"]*self.data.meta_data["FlowEpsilon"] # factor 0.5 left out, see paper by 

class AnalyseQQuartic(FlowAnalyser):
	"""
	Quartic topological charge analysis class.
	"""
	observable_name = r"Topological charge at $Q^4$"
	observable_name_compact = "topq4"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\langle Q^4 \rangle [fm^{-2}]$"

	def __init__(self,*args,**kwargs):
		super(AnalyseQQuartic,self).__init__(*args,**kwargs)
		self.y **= 4

class AnalyseQQ(FlowAnalyser):
	"""
	Quartic topological charge analysis class.
	"""
	observable_name = r"Topological charge at $Q^2$"
	observable_name_compact = "topq4"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\langle Q^2 \rangle [fm^{-1}]$"

	def __init__(self,*args,**kwargs):
		super(AnalyseQQ,self).__init__(*args,**kwargs)
		self.y **= 2


class AnalyseQtQZero(FlowAnalyser):
	"""
	Topological charge QtQ0 analysis class.
	"""
	observable_name = r"Topological charge evolved at flow time $t_0$"
	observable_name_compact = "qtqzero"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\langle Q_{t}Q_{t_0} \rangle [GeV]$"
	observable_output_folder_old = ""

	def __init__(self,*args,**kwargs):
		super(AnalyseQtQZero,self).__init__(*args,**kwargs)
		self.observable_output_folder_path_old = self.observable_output_folder_path
		
	def setQ0(self, q_flow_time_zero, unit_test = False):
		"""
		Sets the flow time we are to analyse for.
		"""
		self.q_flow_time_zero = q_flow_time_zero
		self.flow_time_zero_index = np.argmin(np.abs( self.a * np.sqrt(8*self.x*self.data.meta_data["FlowEpsilon"]) - q_flow_time_zero))

		self.observable_name = "Topological charge evolved at t=%.2f" % (q_flow_time_zero*self.data.meta_data["FlowEpsilon"])

		if unit_test:
			# Performs a deep copy of self.y values(otherwise we will overwrite what we have)
			y_temp1 = copy.deepcopy(self.y)
			y_temp2 = np.zeros(y_temp1.shape)
			y_q0 = copy.deepcopy(self.y[:,self.flow_time_zero_index])

		# Numpy method for matrix-vector multiplication
		self.y = (self.y.T*self.y[:,self.flow_time_zero_index]).T

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old,"t0_%s" % str(self.flow_time_zero_index).strip("."))
		check_folder(self.observable_output_folder_path,self.dryrun,self.verbose)

		if unit_test:
			# Hard-coded matrix-vector multiplication
			for iCfg in xrange(self.N_configurations):
				for iFlow in xrange(self.NFlows):
					y_temp2[iCfg,iFlow] = y_temp1[iCfg,iFlow] * y_q0[iCfg]

			# Matrix-matrix comparison
			for i,j in zip(self.y,y_temp2):
				for ii,jj in zip(i,j):
					if np.abs(ii-jj)>1e-16:
						print "BAD: multiplications do not match."
						exit(1)
			else:
				print "GOOD: multiplications match."

class AnalyseTopologicalSusceptibility(FlowAnalyser):
	"""
	Topological susceptibility analysis class.
	"""
	observable_name = "Topological Susceptibility"
	observable_name_compact = "topsus"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"

	def __init__(self,*args,**kwargs):
		super(AnalyseTopologicalSusceptibility,self).__init__(*args,**kwargs)
		self.y **= 2
		self.set_size()

	def set_size(self): # HARDCODE DIFFERENT SIZES HERE!
		# Retrieves lattice spacing
		self.a = self.getLatticeSpacing(self.data.meta_data["beta"])
		
		# Ugly, hardcoded lattice size setting
		if self.data.meta_data["beta"] == 6.0:
			lattice_size = 24**3*48
			self.chi = ptools._chi_beta6_0
			self.chi_std = ptools._chi_beta6_0_error
			self.chi_derivative = ptools._chi_beta6_0_derivative
		elif self.data.meta_data["beta"] == 6.1:
			lattice_size = 28**3*56
			self.chi = ptools._chi_beta6_1
			self.chi_std = ptools._chi_beta6_1_error
			self.chi_derivative = ptools._chi_beta6_1_derivative
		elif self.data.meta_data["beta"] == 6.2:
			lattice_size = 32**3*64
			self.chi = ptools._chi_beta6_2
			self.chi_std = ptools._chi_beta6_2_error
			self.chi_derivative = ptools._chi_beta6_2_derivative
		elif self.data.meta_data["beta"] == 6.45:
			lattice_size = 48**3*96
			self.chi = ptools._chi_beta6_45
			self.chi_std = ptools._chi_beta6_45_error
			self.chi_derivative = ptools._chi_beta6_45_derivative
		else:
			raise ValueError("Unrecognized beta value: %g" % meta_data["beta"])

		# Sets up constants used in the chi function for topological susceptibility
		self.V = lattice_size
		self.hbarc = 0.19732697 #eV micro m
		self.const = self.hbarc/self.a/self.V**(1./4)

def main(args):
	batch_name = args[0]
	batch_folder = args[1]
	DirectoryList = GetDirectoryTree(batch_name,batch_folder)
	N_bs = 500
	dryrun = False
	verbose = True
	parallel = True
	numprocs = 8
	create_perflow_data = False
	# print DirectoryList

	# Analysis timers
	pre_time = time.clock()
	observable_strings = []

	# Analyses plaquette data if present in arguments
	if 'plaq' in args:
		plaq_analysis = AnalysePlaquette(DirectoryList, batch_name, dryrun = dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose, create_perflow_data = create_perflow_data)
		plaq_analysis.boot(N_bs)
		plaq_analysis.jackknife()
		# plaq_analysis.y_limits = [0.55,1.05]
		plaq_analysis.plot_original()
		plaq_analysis.plot_boot()
		plaq_analysis.plot_jackknife()
		plaq_analysis.autocorrelation()
		plaq_analysis.plot_autocorrelation(0)
		plaq_analysis.plot_autocorrelation(-1)
		plaq_analysis.plot_mc_history(0)
		plaq_analysis.plot_mc_history(-1)
		plaq_analysis.plot_original()
		plaq_analysis.plot_boot()
		plaq_analysis.plot_jackknife()
		plaq_analysis.plot_histogram(0,r"$P_{\mu\nu}$",x_limits='auto')
		plaq_analysis.plot_histogram(-1,r"$P_{\mu\nu}$",x_limits='auto')
		plaq_analysis.plot_integrated_correlation_time()
		plaq_analysis.plot_integrated_correlation_time()

		if plaq_analysis.bootstrap_performed and plaq_analysis.autocorrelation_performed:
			write_data_to_file(plaq_analysis,dryrun = dryrun)

	if 'topc' in args or 'topsus' in args or 'topcq4' in args or 'qtqzero' in args:
		topc_analysis = AnalyseTopologicalCharge(DirectoryList, batch_name, dryrun = dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose, create_perflow_data = create_perflow_data)
		if 'topc' in args:

			topc_analysis.boot(N_bs)
			topc_analysis.jackknife()
			# topc_analysis.y_limits = [-10,10]
			topc_analysis.plot_original()
			topc_analysis.plot_boot()
			topc_analysis.plot_jackknife()
			topc_analysis.autocorrelation()
			topc_analysis.plot_autocorrelation(0)
			topc_analysis.plot_autocorrelation(-1)

			topc_analysis.plot_mc_history(0)
			topc_analysis.plot_mc_history(-1)

			# topc_analysis.plot_boot()
			# topc_analysis.plot_original()
			# topc_analysis.plot_jackknife()
			# topc_analysis.plot_histogram(0,r"$Q$[GeV]",x_limits='auto')
			# topc_analysis.plot_histogram(-1,r"$Q$[GeV]",x_limits='auto')
			# topc_analysis.plot_integrated_correlation_time()
			# topc_analysis.plot_integrated_correlation_time()

			if topc_analysis.bootstrap_performed and topc_analysis.autocorrelation_performed:
				write_data_to_file(topc_analysis,dryrun = dryrun)

		if 'topcq4' in args:
			topcq4_analysis = AnalyseQQuartic(DirectoryList, batch_name, data = topc_analysis.data, dryrun = dryrun, parallel = parallel, numprocs = numprocs, verbose=verbose)
			topcq4_analysis.boot(N_bs)
			topcq4_analysis.jackknife()
			topcq4_analysis.plot_original()
			topcq4_analysis.plot_jackknife()
			topcq4_analysis.plot_boot()
			topcq4_analysis.autocorrelation()
			topcq4_analysis.plot_autocorrelation(-1)
			topcq4_analysis.plot_autocorrelation(0)
			topcq4_analysis.plot_mc_history(0)
			topcq4_analysis.plot_mc_history(-1)
			topcq4_analysis.plot_original()
			topcq4_analysis.plot_boot()
			topcq4_analysis.plot_jackknife()
			topcq4_analysis.plot_histogram(0,topcq4_analysis.y_label)
			topcq4_analysis.plot_histogram(-1,topcq4_analysis.y_label)
			topcq4_analysis.plot_integrated_correlation_time()
			topcq4_analysis.plot_integrated_correlation_time()

		if 'topcqq' in args:
			topcqq_analysis = AnalyseQQ(DirectoryList, batch_name, data = topc_analysis.data, dryrun = dryrun, parallel = parallel, numprocs = numprocs, verbose=verbose)
			topcqq_analysis.boot(N_bs)
			topcqq_analysis.jackknife()
			topcqq_analysis.plot_original()
			topcqq_analysis.plot_jackknife()
			topcqq_analysis.plot_boot()
			topcqq_analysis.autocorrelation()
			topcqq_analysis.plot_autocorrelation(-1)
			topcqq_analysis.plot_autocorrelation(0)
			topcqq_analysis.plot_mc_history(0)
			topcqq_analysis.plot_mc_history(-1)
			topcqq_analysis.plot_original()
			topcqq_analysis.plot_boot()
			topcqq_analysis.plot_jackknife()
			topcqq_analysis.plot_histogram(0,topcqq_analysis.y_label)
			topcqq_analysis.plot_histogram(-1,topcqq_analysis.y_label)
			topcqq_analysis.plot_integrated_correlation_time()
			topcqq_analysis.plot_integrated_correlation_time()			

		if 'qtqzero' in args:
			qzero_flow_times = [0.1,0.2,0.3,0.4,0.5,0.6]
			qtqzero_analysis = AnalyseQtQZero(DirectoryList, batch_name, data = topc_analysis.data, dryrun = dryrun, parallel = parallel, numprocs = numprocs, verbose=verbose)
			for qzero_flow_time in qzero_flow_times:
				qtqzero_analysis.setQ0(qzero_flow_time)
				qtqzero_analysis.boot(N_bs)
				qtqzero_analysis.jackknife()
				qtqzero_analysis.plot_original()
				qtqzero_analysis.plot_jackknife()
				qtqzero_analysis.plot_boot()

		if 'topsus' in args:
			topsus_analysis = AnalyseTopologicalSusceptibility(DirectoryList, batch_name, dryrun = dryrun, data=topc_analysis.data, parallel=parallel, numprocs=numprocs, verbose=verbose)
			topsus_analysis.boot(N_bs,F = topsus_analysis.chi, F_error = topsus_analysis.chi_std, store_raw_bs_values = True)
			topsus_analysis.jackknife(F = topsus_analysis.chi, F_error = topsus_analysis.chi_std, store_raw_jk_values = True)
			# topsus_analysis.y_limits  = [0.05,0.5]
			topsus_analysis.plot_original()
			topsus_analysis.plot_boot()
			topsus_analysis.plot_jackknife()
			topsus_analysis.autocorrelation()
			topsus_analysis.plot_autocorrelation(0)
			topsus_analysis.plot_autocorrelation(-1)
			topsus_analysis.plot_mc_history(0)
			topsus_analysis.plot_mc_history(-1)
			topsus_analysis.plot_original()
			topsus_analysis.plot_boot()
			topsus_analysis.plot_jackknife()
			topsus_analysis.plot_histogram(0,r"$\chi^{1/4}$[GeV]",x_limits='auto')
			topsus_analysis.plot_histogram(-1,r"$\chi^{1/4}$[GeV]",x_limits='auto')
			topsus_analysis.plot_integrated_correlation_time()
			topsus_analysis.plot_integrated_correlation_time()

			if topsus_analysis.bootstrap_performed and topsus_analysis.autocorrelation_performed:
				write_data_to_file(topsus_analysis, dryrun = dryrun)

	if 'energy' in args:
		r0 = 0.5
		energy_analysis = AnalyseEnergy(DirectoryList, batch_name, dryrun = dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose, create_perflow_data = create_perflow_data)
		energy_analysis.boot(N_bs,store_raw_bs_values=True)
		energy_analysis.jackknife(store_raw_jk_values=True)
		x_values = energy_analysis.data.meta_data["FlowEpsilon"] * energy_analysis.x / r0**2 * energy_analysis.a**2
		# energy_analysis.y_limits = [ 0 , np.max(energy_analysis.correction_function(energy_analysis.unanalyzed_y))]
		energy_analysis.plot_original(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_boot(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_jackknife(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.autocorrelation()
		energy_analysis.plot_autocorrelation(0)
		energy_analysis.plot_autocorrelation(-1)
		energy_analysis.plot_mc_history(0)
		energy_analysis.plot_mc_history(-1)
		energy_analysis.plot_original(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_boot(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_jackknife(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_histogram(0,r"$E$[GeV]",x_limits='auto')
		energy_analysis.plot_histogram(-1,r"$E$[GeV]",x_limits='auto')
		energy_analysis.plot_integrated_correlation_time()
		energy_analysis.plot_integrated_correlation_time()

		if energy_analysis.bootstrap_performed and energy_analysis.autocorrelation_performed:
			write_data_to_file(energy_analysis, dryrun = dryrun)

	post_time = time.clock()
	print "="*100
	print "Analysis of batch %s observables %s in %.2f seconds" % (batch_name, ", ".join([i.lower() for i in args[2:]]), (post_time-pre_time))
	print "="*100

if __name__ == '__main__':
	if not sys.argv[1:]:
		# args = [['prodRunBeta6_0','output','plaq','topc','energy','topsus'],
		# 		['prodRunBeta6_1','output','plaq','topc','energy','topsus']]

		# args = [['beta6_0','data','plaq','topc','energy','topsus'],
		# 		['beta6_1','data','plaq','topc','energy','topsus'],
		# 		['beta6_2','data','plaq','topc','energy','topsus']]

		# args = [['beta6_0','data2','plaq','topc','energy','topsus'],
		# 		['beta6_1','data2','plaq','topc','energy','topsus'],
		# 		['beta6_2','data2','plaq','topc','energy','topsus']]

		# args = [['beta6_0','data4','plaq','topc','energy','topsus'],
		# 		['beta6_1','data4','plaq','topc','energy','topsus'],
		# 		['beta6_2','data4','plaq','topc','energy','topsus']]

		# args = [['beta6_0','data2','plaq','topc','energy','topsus','qtqzero','topcq4','topcqq'],
		# 		['beta6_1','data2','plaq','topc','energy','topsus','qtqzero','topcq4','topcqq'],
		# 		['beta6_2','data2','plaq','topc','energy','topsus','qtqzero','topcq4','topcqq']]

		args = [['beta6_0','data4','topc']]

		# args = [['beta6_2','data3','topcq4']]
		# args = [['beta6_2','data3','qtqzero']]

		# args = [['beta6_1','data','topsus'],
		# 		['beta6_1','data2','topsus'],
		# 		['beta6_1','data3','topsus']]

		# args = [['beta6_0','data2','topsus','energy'],
		# 		['beta6_1','data2','topsus','energy'],
		# 		['beta6_2','data2','topsus','energy']]

		# args = [['beta6_1','data2','topsus']]

		# args = [['beta6_0','data','topc','topsus','energy']]

		# args = [['test_run_new_counting','output','topc','plaq','energy','topsus']]

		# args = [['beta6_1','data','plaq','topc','energy','topsus']]
		# args = [['beta6_2','data','energy']]

		for a in args:
			main(a)
	else:
		main(sys.argv[1:])

	# plt.show()