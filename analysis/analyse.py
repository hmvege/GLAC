from tools.folderreadingtools import GetDirectoryTree, GetFolderContents, write_data_to_file
from statistics.jackknife import Jackknife
from statistics.bootstrap import Bootstrap
from statistics.autocorrelation import Autocorrelation
import os, numpy as np, matplotlib.pyplot as plt, sys, pandas as pd, multiprocessing, pickle, time

# from tqdm import trange

#### Parallel helper functions ####
def _autocorrelation_parallel_core(input_values):
	ac = Autocorrelation(input_values[0],use_numpy=input_values[1],F=input_values[2])
	return ac.R, ac.R_error, ac.integrated_autocorrelation_time(), ac.integrated_autocorrelation_time_error()

def _bootstrap_parallel_core(input_values):
	data, N_bs, bs_statistic, F, non_bs_stats, index_lists = input_values
	bs = Bootstrap(data, N_bs, bootstrap_statistics = bs_statistic, F = F, non_bs_stats = non_bs_stats, index_lists = index_lists)
	return bs.bs_avg, bs.bs_std, bs.avg_original, bs.std_original, bs.bs_data, bs.data_original

def _jackknife_parallel_core(input_values):
	data, F, jk_statistics, non_jk_statistics = input_values
	jk = Jackknife(data, F = F, jk_statistics = jk_statistics, non_jk_statistics = non_jk_statistics)
	return jk.jk_avg, jk.jk_std, jk.jk_data

def _default_return(x):
	# For use instead of lambda x : x in parallel
	return x

def _default_return_squared(x):
	# For use instead of lambda x**2 : x**2 in parallel
	return x*x

class FlowAnalyser(object):
	observable_name = "Missing_Observable_Name"
	observable_name_compact = "missing_obs_name"
	x_label = "Missing x-label"
	y_label = "Missing y-label"
	mark_interval = 5
	error_mark_interval = 5
	autocorrelations_limits = 1
	dpi = None

	def __init__(self,files,batch_name,data=None,dryrun=False,flow=True, parallel = False, numprocs = 4):
		# Sets up global constants
		self.batch_name = batch_name
		self.N_bs = None
		self.dryrun = dryrun
		self.flow = flow

		# Parallel variables
		self.parallel = parallel
		self.numprocs = numprocs

		print "="*100
		print "Batch:       %s" % self.batch_name
		print "Observables: %s" % self.observable_name
		print "="*100

		# Enables possibility of providing data
		if data == None:
			self.data = GetFolderContents(files,flow=self.flow)
		else:
			self.data = data

		# print "Size of data: %.4g kB" % (sys.getsizeof(self.data.data_y)/1024.0)

		# Sets up variables
		self.y = self.data.data_y
		self.x = self.data.data_x

		# Max plotting window variables
		self.y_limits = []

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

	def boot(self,N_bs,bs_statistic = np.mean, F = _default_return, non_bs_stats = _default_return):
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
								[bs_statistic for i in xrange(self.NFlows)],
								[F for i in xrange(self.NFlows)],
								[non_bs_stats for i in xrange(self.NFlows)],
								[index_lists for i in xrange(self.NFlows)])

			# Initializes multiprocessing
			pool = multiprocessing.Pool(processes=self.numprocs)								

			# Runs parallel processes. Can this be done more efficiently?
			results = pool.map(_bootstrap_parallel_core,input_values)

			# Garbage collection for multiprocessing instance
			pool.close()

			# Populating bootstrap data
			for i in xrange(self.NFlows):
				self.bs_y[i] = results[i][0]
				self.bs_y_std[i] = results[i][1]
				self.unanalyzed_y[i] = results[i][2]
				self.unanalyzed_y_std[i] = results[i][3]

				# Stores last data for plotting in histogram later
				self.bs_y_data[i] = results[i][4]
				self.unanalyzed_y_data[i] = results[i][5]

		else:
			print "Not running parallel bootstrap for %s!" % self.observable_name
			# Non-parallel method for calculating jackknife
			for i in xrange(self.NFlows):
				bs = Bootstrap(self.y[:,i], N_bs, bootstrap_statistics = bs_statistic, F = F, non_bs_stats = non_bs_stats, index_lists = index_lists)
				self.bs_y[i] = bs.bs_avg
				self.bs_y_std[i] = bs.bs_std
				self.unanalyzed_y[i] = bs.avg_original
				self.unanalyzed_y_std[i] = bs.std_original

				# Stores last data for plotting in histogram later
				self.bs_y_data[i] = bs.bs_data
				self.unanalyzed_y_data[i] = bs.data_original

		# Sets performed flag to true
		self.bootstrap_performed = True

	def jackknife(self, jk_statistics = np.average, F = _default_return, non_jk_statistics = _default_return):
		if self.parallel:
			# Sets up jobs for parallel processing
			input_values = zip(	[self.y[:,i] for i in xrange(self.NFlows)],
								[F for i in xrange(self.NFlows)],
								[jk_statistics for i in xrange(self.NFlows)],
								[non_jk_statistics for i in xrange(self.NFlows)])

			# Prints job size
			# print sys.getsizeof(input_values)/1024.0, "kB"

			# Initializes multiprocessing
			pool = multiprocessing.Pool(processes=self.numprocs)								

			# Runs parallel processes. Can this be done more efficiently?
			results = pool.map(_jackknife_parallel_core,input_values)

			# Closes multiprocessing instance for garbage collection
			pool.close()

			# Populating jackknife results
			for i in xrange(self.NFlows):
				self.jk_y[i] = results[i][0]
				self.jk_y_std[i] = results[i][1]
				self.jk_y_data[i] = results[i][2]
		else:
			print "Not running parallel jackknife for %s!" % self.observable_name
			# Non-parallel method for calculating jackknife
			for i in xrange(self.NFlows):
				jk = Jackknife(self.y[:,i],F = F, jk_statistics = jk_statistics, non_jk_statistics = non_jk_statistics)
				self.jk_y[i] = jk.jk_avg
				self.jk_y_std[i] = jk.jk_std
				self.jk_y_data[i] = jk.jk_data

		# Sets performed flag to true
		self.jackknife_performed = True

	def autocorrelation(self,use_numpy=True, F = _default_return):
		# Gets autocorrelation
		if self.parallel:
			# Sets up jobs for parallel processing
			input_values = zip(	[self.y[:,i] for i in xrange(self.NFlows)],
								[use_numpy for i in xrange(self.NFlows)],
								[F for i in xrange(self.NFlows)])

			# Sets up parallel job
			pool = multiprocessing.Pool(processes=self.numprocs)

			# Initiates parallel jobs
			results = pool.map(_autocorrelation_parallel_core,input_values)

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
			print "Not running parallel autocorrelation for %s!" % self.observable_name
			# Non-parallel method for calculating autocorrelation
			for i in xrange(self.NFlows):
				ac = Autocorrelation(self.y[:,i],use_numpy=use_numpy,F=F)
				self.autocorrelations[i] = ac.R
				self.autocorrelations_errors[i] = ac.R_error
				self.integrated_autocorrelation_time[i] = ac.integrated_autocorrelation_time()
				self.integrated_autocorrelation_time_error[i] = ac.integrated_autocorrelation_time_error()
				self.autocorrelation_error_correction[i] = np.sqrt(2*self.integrated_autocorrelation_time[i])

				# Small progressbar
				sys.stdout.write("\rCalculating autocorrelation: %4.1f%% done" % (100*float(i)/float(self.NFlows)))
				sys.stdout.flush()

			# Finalizes
			sys.stdout.write("\rCalculating autocorrelation: 100.0%% done")
			sys.stdout.flush()

		# print "PRINTING IN Autocorrelation @ analyse.py"
		# print self.integrated_autocorrelation_time
		# print self.integrated_autocorrelation_time_error


		# Sets performed flag to true
		self.autocorrelation_performed = True

	def plot_jackknife(self, x = None, correction_function = lambda x : x, fname_addon = ""):
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
		fname = "../figures/{0:<s}/{1:<s}_jackknife_beta{2:<s}{3:<s}.png".format(self.batch_name,self.observable_name_compact,str(self.beta).replace('.','_'),fname_addon)

		# Plots the jackknifed data
		self.__plot_error_core(x,correction_function(y),correction_function(y_std),title_string,fname)

	def plot_boot(self,plot_bs=True, x = None, correction_function = lambda x : x, fname_addon = ""):
		# Checks that the bootstrap has been performed.
		if not self.bootstrap_performed and plot_bs:
			raise ValueError("Bootstrap has not been performed yet.")

		# Retrieves relevant data and sets up the arrays to be plotted
		if type(x) != np.ndarray:
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x*self.data.meta_data["FlowEpsilon"])
		if plot_bs:
			y = self.bs_y
			y_std = self.bs_y_std*self.autocorrelation_error_correction
		else:
			y = self.unanalyzed_y
			y_std = self.unanalyzed_y_std*self.autocorrelation_error_correction
			# print y_std 

		# Sets up the title and filename strings
		if plot_bs:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta,self.N_bs)
			fname = "../figures/{0:<s}/{1:<s}_bootstrap_Nbs{3:<d}_beta{2:<s}{4:<s}.png".format(self.batch_name,self.observable_name_compact,str(self.beta).replace('.','_'),self.N_bs,fname_addon)
		else:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name, self.data.meta_data["NFlows"],self.beta)
			fname = "../figures/{0:<s}/{1:<s}_original_beta{2:<s}{3:<s}.png".format(self.batch_name,self.observable_name_compact,str(self.beta).replace('.','_'),fname_addon)

		# Plots either bootstrapped or regular stuff
		self.__plot_error_core(x,correction_function(y),correction_function(y_std),title_string,fname)

	def plot_original(self, x = None, correction_function = lambda x : x, fname_addon = ""):
		"""
		Plots the default analysis, mean and std of the observable.
		"""
		self.plot_boot(plot_bs=False, x = x, correction_function = correction_function, fname_addon = fname_addon)

	def __plot_error_core(self,x,y,y_std,title_string,fname):
		# Plots the jackknifed data
		fig = plt.figure()
		ax = fig.add_subplot(111)

		ax.errorbar(x,y,yerr=y_std,fmt=".",color="0",ecolor="r",label=self.observable_name,markevery=self.mark_interval,errorevery=self.error_mark_interval)
		# pl.plot(x, y, 'k', color='#CC4F1B')
		# pl.fill_between(x, y-error, y+error,
		#     alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

		if len(self.y_limits) == 0:
			self.y_limits = [-np.max(y),np.max(y)]

		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)
		ax.set_ylim(self.y_limits)
		ax.grid(True)
		ax.set_title(title_string)
		if not self.dryrun: 
			fig.savefig(fname,dpi=self.dpi)
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
		# for i in xrange(N_autocorr):
			# y[i] = np.max(self.autocorrelations[:,i])
		if plot_abs_value:
			y = np.abs(self.autocorrelations[flow_time,:])
		else:
			y = self.autocorrelations[flow_time,:]
		
		y_std = self.autocorrelations_errors[flow_time,:]

		# Sets up the title and filename strings
		title_string = r"Autocorrelation of %s at flow time $t=%.2f$, $\beta=%.2f$, $N_{cfg}=%2d$" % (self.observable_name,flow_time*self.data.meta_data["FlowEpsilon"],self.beta,self.N_configurations)
		fname = "../figures/{0:<s}/{1:<s}_autocorrelation_flowt{2:<d}_beta{3:<s}.png".format(self.batch_name,self.observable_name_compact,flow_time, str(self.beta).replace('.','_'))

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
			fig.savefig(fname,dpi=self.dpi)
		print "Figure created in %s" % fname
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
		fname = "../figures/{0:<s}/{1:<s}_histogram_flowt{2:<d}_beta{3:<s}.png".format(self.batch_name,self.observable_name_compact,abs(flow_time),str(self.beta).replace('.','_'))

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
			plt.savefig(fname)
		print "Figure created in %s" % fname

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
		fname = "../figures/{0:<s}/{1:<s}_mchistory_flowt{2:<d}_beta{3:<s}.png".format(self.batch_name,self.observable_name_compact,flow_time,str(self.beta).replace('.','_'))

		# Sets up plot
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(correction_function(self.unanalyzed_y_data[flow_time]),color="0",label=self.observable_name)
		ax.set_xlabel(r"Monte Carlo time")
		ax.set_ylabel(r"")
		ax.set_title(title_string)
		ax.grid(True)
		ax.legend()
		if not self.dryrun: 
			fig.savefig(fname,dpi=self.dpi)
		print "Figure created in %s" % fname
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
	observable_name = "Topological Charge"
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

	def correction_function(self, y):
		return -y*self.x*self.x*self.data.meta_data["FlowEpsilon"]*self.data.meta_data["FlowEpsilon"]/64.0 # factor 0.5 left out, see paper by 

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
		self.set_size()

	def set_size(self):
		self.a = self.getLatticeSpacing(self.data.meta_data["beta"])
		
		# Ugly, hardcoded lattice size setting
		if self.data.meta_data["beta"] == 6.0:
			lattice_size = 24**3*48
		elif self.data.meta_data["beta"] == 6.1:
			lattice_size = 28**3*56
		elif self.data.meta_data["beta"] == 6.2:
			lattice_size = 32**3*64
		elif self.data.meta_data["beta"] == 6.45:
			lattice_size = 48**3*96
		else:
			raise ValueError("Unrecognized beta value: %g" % meta_data["beta"])

		# Sets up constants used in the chi function for topological susceptibility
		self.V = lattice_size
		self.hbarc = 0.19732697 #eV micro m
		self.const = self.hbarc/self.a/self.V**(1./4)

	def jackknife(self,*args,**kwargs):
		_temp_parallel = self.parallel
		self.parallel = False
		super(AnalyseTopologicalSusceptibility,self).jackknife(*args,**kwargs)
		self.parallel = _temp_parallel

	def boot(self,*args,**kwargs):
		_temp_parallel = self.parallel
		self.parallel = False
		super(AnalyseTopologicalSusceptibility,self).boot(*args,**kwargs)
		self.parallel = _temp_parallel

	def chi(self, Q_squared):
		# Q should be averaged
		return self.const*Q_squared**(1./4)

	@staticmethod
	def stat(x,axis=None):
		return np.mean(x**2,axis=axis)

	@staticmethod
	def return_x_squared(x):
		return x**2

def main(args):
	DirectoryList = GetDirectoryTree(args[0],output_folder=args[1])
	N_bs = 500
	dryrun = False
	parallel = True
	numprocs = 8
	use_numpy_in_autocorrelation = True
	# print DirectoryList

	# Analysis timers
	pre_time = time.clock()
	observable_strings = []

	# Analyses plaquette data if present in arguments
	if 'plaq' in args:
		plaq_analysis = AnalysePlaquette(DirectoryList.getFlow("plaq"), args[0], dryrun = dryrun, parallel=parallel, numprocs=numprocs)
		plaq_analysis.boot(N_bs)
		plaq_analysis.jackknife()
		plaq_analysis.y_limits = [0.55,1.05]
		plaq_analysis.plot_original(fname_addon = "_noErrorCorrection")
		plaq_analysis.plot_boot(fname_addon = "_noErrorCorrection")
		plaq_analysis.plot_jackknife(fname_addon = "_noErrorCorrection")
		plaq_analysis.autocorrelation(use_numpy=use_numpy_in_autocorrelation)
		plaq_analysis.plot_autocorrelation(0)
		plaq_analysis.plot_autocorrelation(-1)
		plaq_analysis.plot_mc_history(0)
		plaq_analysis.plot_mc_history(-1)
		plaq_analysis.plot_original()
		plaq_analysis.plot_boot()
		plaq_analysis.plot_jackknife()
		plaq_analysis.plot_histogram(0,r"$P_{\mu\nu}$",x_limits='analysis')
		plaq_analysis.plot_histogram(-1,r"$P_{\mu\nu}$",x_limits='analysis')

		if plaq_analysis.bootstrap_performed and plaq_analysis.autocorrelation_performed:
			write_data_to_file(plaq_analysis,dryrun = dryrun)

	if 'topc' in args or 'topsus' in args:
		topc_analysis = AnalyseTopologicalCharge(DirectoryList.getFlow("topc"), args[0], dryrun = dryrun, parallel=parallel, numprocs=numprocs)
		if 'topc' in args:
			topc_analysis.boot(N_bs)
			topc_analysis.jackknife()
			topc_analysis.y_limits = [-10,10]
			topc_analysis.plot_original(fname_addon = "_noErrorCorrection")
			topc_analysis.plot_boot(fname_addon = "_noErrorCorrection")
			topc_analysis.plot_jackknife(fname_addon = "_noErrorCorrection")
			topc_analysis.autocorrelation(use_numpy=use_numpy_in_autocorrelation)
			topc_analysis.plot_autocorrelation(0)
			topc_analysis.plot_autocorrelation(-1)
			topc_analysis.plot_mc_history(0)
			topc_analysis.plot_mc_history(-1)
			topc_analysis.plot_boot()
			topc_analysis.plot_original()
			topc_analysis.plot_jackknife()
			topc_analysis.plot_histogram(0,r"$Q$[GeV]",x_limits='analysis')
			topc_analysis.plot_histogram(-1,r"$Q$[GeV]",x_limits='analysis')

			if topc_analysis.bootstrap_performed and topc_analysis.autocorrelation_performed:
				write_data_to_file(topc_analysis,dryrun = dryrun)

		if 'topsus' in args:
			topsus_analysis = AnalyseTopologicalSusceptibility(DirectoryList.getFlow("topc"), args[0], dryrun = dryrun, data=topc_analysis.data, parallel=parallel, numprocs=numprocs)
			topsus_analysis.boot(N_bs,bs_statistic = topsus_analysis.stat, F = topsus_analysis.chi, non_bs_stats = topsus_analysis.return_x_squared)
			topsus_analysis.jackknife(jk_statistics = topsus_analysis.stat, F = topsus_analysis.chi, non_jk_statistics = topsus_analysis.return_x_squared)
			topsus_analysis.y_limits = [0.05,0.24]
			topsus_analysis.plot_original(fname_addon = "_noErrorCorrection")
			topsus_analysis.plot_boot(fname_addon = "_noErrorCorrection")
			topsus_analysis.plot_jackknife(fname_addon = "_noErrorCorrection")
			topsus_analysis.autocorrelation(use_numpy=use_numpy_in_autocorrelation, F = _default_return_squared) # Dosen't make sense to do the autocorrelation of the topoligical susceptibility since it is based on a data mean
			topsus_analysis.plot_autocorrelation(0)
			topsus_analysis.plot_autocorrelation(-1)
			topsus_analysis.plot_mc_history(0)
			topsus_analysis.plot_mc_history(-1)
			topsus_analysis.plot_original()
			topsus_analysis.plot_boot()
			topsus_analysis.plot_jackknife()
			topsus_analysis.plot_histogram(0,r"$\chi^{1/4}$[GeV]",x_limits='analysis')
			topsus_analysis.plot_histogram(-1,r"$\chi^{1/4}$[GeV]",x_limits='analysis')

			if topsus_analysis.bootstrap_performed and topsus_analysis.autocorrelation_performed:
				write_data_to_file(topsus_analysis,dryrun = dryrun)

	if 'energy' in args:
		r0 = 0.5
		energy_analysis = AnalyseEnergy(DirectoryList.getFlow("energy"), args[0], dryrun = dryrun, parallel=parallel, numprocs=numprocs)
		energy_analysis.boot(N_bs)
		energy_analysis.jackknife()
		x_values = energy_analysis.data.meta_data["FlowEpsilon"] * energy_analysis.x / r0**2 * energy_analysis.a**2
		energy_analysis.y_limits = [ 0 , np.max(energy_analysis.correction_function(energy_analysis.unanalyzed_y))]
		energy_analysis.plot_original(x = x_values, correction_function = energy_analysis.correction_function, fname_addon = "_noErrorCorrection")
		energy_analysis.plot_boot(x = x_values, correction_function = energy_analysis.correction_function, fname_addon = "_noErrorCorrection")
		energy_analysis.plot_jackknife(x = x_values, correction_function = energy_analysis.correction_function, fname_addon = "_noErrorCorrection")
		energy_analysis.autocorrelation(use_numpy=use_numpy_in_autocorrelation)
		energy_analysis.plot_autocorrelation(0)
		energy_analysis.plot_autocorrelation(-1)
		energy_analysis.plot_mc_history(0,correction_function = lambda x : -x/64.0)
		energy_analysis.plot_mc_history(-1,correction_function = lambda x : -x/64.0)
		energy_analysis.plot_original(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_boot(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_jackknife(x = x_values, correction_function = energy_analysis.correction_function)
		energy_analysis.plot_histogram(0,r"$E$[GeV]",x_limits='analysis')
		energy_analysis.plot_histogram(-1,r"$E$[GeV]",x_limits='analysis')

		if energy_analysis.bootstrap_performed and energy_analysis.autocorrelation_performed:
			write_data_to_file(energy_analysis,dryrun = dryrun)

	post_time = time.clock()
	print "="*100
	print "Analysis of batch %s observables %s in %.2f seconds" % (args[0], ", ".join([i.lower() for i in args[1:]]), (post_time-pre_time))
	print "="*100

if __name__ == '__main__':
	if not sys.argv[1:]:
		# args = [['prodRunBeta6_0','output','plaq','topc','energy','topsus'],
		# 		['prodRunBeta6_1','output','plaq','topc','energy','topsus']]

		# args = [['beta6_0','data','plaq','topc','energy','topsus'],
		# 		['beta6_1','data','plaq','topc','energy','topsus'],
		# 		['beta6_2','data','plaq','topc','energy','topsus']]

		args = [['beta6_0','data2','plaq','topc','energy','topsus'],
				['beta6_1','data2','plaq','topc','energy','topsus'],
				['beta6_2','data2','plaq','topc','energy','topsus']]

		# args = [['beta6_2','data2','plaq','topc','energy','topsus']]

		# args = [['beta6_0','data2','topsus'],
		# 		['beta6_1','data2','topsus'],
		# 		['beta6_2','data2','topsus']]

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