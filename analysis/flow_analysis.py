from flow_analyser import FlowAnalyser
from tools.folderreadingtools import DataReader
import statistics.parallel_tools as ptools
import os
import numpy as np
import copy
import sys
import time

# from tqdm import trange

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
	scaling_factor = 1.0 # In case the data set contains unscaled x values
	x_label = r"$t/r_0^2$" # Dimensionsless, Implied multiplication by a^2
	y_label = r"$t^2\langle E \rangle$" # Energy is dimension 4, while t^2 is dimension invsere 4, or length/time which is inverse energy, see Peskin and Schroeder

	def __init__(self, *args, **kwargs):
		super(AnalyseEnergy, self).__init__(*args, **kwargs)
		self.y *= -1.0/64.0

	def correction_function(self, y):
		# *  self.flow_epsilon * self.flow_epsilon
		return y * self.x * self.x * self.scaling_factor # factor 0.5 left out, see paper by 

class AnalyseQQuartic(FlowAnalyser):
	"""
	Quartic topological charge analysis class.
	"""
	observable_name = r"Topological charge at $Q^4$"
	observable_name_compact = "topq4"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\langle Q^4 \rangle [fm^{-2}]$"

	def __init__(self, *args, **kwargs):
		super(AnalyseQQuartic, self).__init__(*args, **kwargs)
		self.y **= 4

class AnalyseQtQZero(FlowAnalyser):
	"""
	Topological charge QtQ0 analysis class.
	"""
	observable_name = r"Topological charge evolved at flow time $t_0$"
	observable_name_compact = "qtqzero"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\langle Q_{t}Q_{t_0} \rangle [GeV]$" # $\chi_t^{1/4}[GeV]$
	observable_output_folder_old = ""

	def __init__(self, *args, **kwargs):
		super(AnalyseQtQZero, self).__init__(*args, **kwargs)
		self.observable_output_folder_path_old = self.observable_output_folder_path
		self.y_original = copy.deepcopy(self.y)

	def setQ0(self, q_flow_time_zero_percent, unit_test=False, y_label=None):
		"""
		Sets the flow time we are to analyse for
		q_flow_time_zero_percent: float between 0.0 and 1.0, in which we choose what percentage point of the data we set as q0.
		E.g. if it is 0.9, it will be the Q that is closest to 90% of the whole flowed time
		"""

		# self.set_size()

		# Finds the q flow time zero value
		self.q_flow_time_zero = q_flow_time_zero_percent * (self.a * np.sqrt(8*self.x))[-1]
		
		# Finds the flow time zero index
		self.flow_time_zero_index = np.argmin(np.abs(self.a * np.sqrt(8*self.x) - self.q_flow_time_zero))
		
		# Sets file name
		self.observable_name = "Topological charge evolved at t=%.3f" % (self.q_flow_time_zero)

		if unit_test:
			# Performs a deep copy of self.y values(otherwise we will overwrite what we have)
			y_temp1 = copy.deepcopy(self.y_original)
			y_temp2 = np.zeros(y_temp1.shape)

		# Manual method for multiplying the matrices
		y_q0 = copy.deepcopy(self.y_original[:,self.flow_time_zero_index])
		self.y = copy.deepcopy(self.y_original)

		# Multiplying QtQ0
		for iFlow in xrange(self.y.shape[1]):
			self.y[:,iFlow] *= y_q0

		if y_label != None:
			self.y_label = y_label

		# Creates a new folder to store t0 results in
		self.observable_output_folder_path = os.path.join(self.observable_output_folder_path_old, "t0_%s" % str(self.flow_time_zero_index).strip("."))
		check_folder(self.observable_output_folder_path, self.dryrun, self.verbose)

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


	# def set_size(self): # HARDCODE DIFFERENT SIZES HERE!
	# 	# Retrieves lattice spacing
	# 	self.a = self.getLatticeSpacing(self.beta)
		
	# 	# Ugly, hardcoded lattice size setting
	# 	if self.beta == 6.0:
	# 		lattice_size = 24**3*48
	# 		self.chi = ptools._chi_beta6_0
	# 		self.chi_std = ptools._chi_beta6_0_error
	# 		self.function_derivative = ptools._chi_beta6_0_derivative
	# 	elif self.beta == 6.1:
	# 		lattice_size = 28**3*56
	# 		self.chi = ptools._chi_beta6_1
	# 		self.chi_std = ptools._chi_beta6_1_error
	# 		self.function_derivative = ptools._chi_beta6_1_derivative
	# 	elif self.beta == 6.2:
	# 		lattice_size = 32**3*64
	# 		self.chi = ptools._chi_beta6_2
	# 		self.chi_std = ptools._chi_beta6_2_error
	# 		self.function_derivative = ptools._chi_beta6_2_derivative
	# 	elif self.beta == 6.45:
	# 		lattice_size = 48**3*96
	# 		self.chi = ptools._chi_beta6_45
	# 		self.chi_std = ptools._chi_beta6_45_error
	# 		self.function_derivative = ptools._chi_beta6_45_derivative
	# 	else:
	# 		raise ValueError("Unrecognized beta value: %g" % meta_data["beta"])

	# 	# Sets up constants used in the chi function for topological susceptibility
	# 	self.V = lattice_size
	# 	self.hbarc = 0.19732697 #eV micro m
	# 	self.const = self.hbarc/self.a/self.V**(1./4)


class AnalyseTopologicalChargeInEuclideanTime(FlowAnalyser):
	"""
	Analysis of the topological charge in Euclidean Time.
	"""
	# overwrite_meta_data_flows = True
	observable_name = "Topological Charge in Euclidean Time"
	observable_name_compact = "topct"
	x_label = r"$t_{Euclidean}[fm]$"
	y_label = r"$\langle Q^2 \rangle_{Euclidean}[GeV]$"

	def __init__(self, *args, **kwargs):
		super(AnalyseTopologicalChargeInEuclideanTime, self).__init__(*args, **kwargs)
		self.y = (self.y**2)**(0.25)


		print "MAKE AnalyseTopologicalChargeInEuclideanTime similar to the QtQ0!! --> exiting"
		exit(1) 

def main(args):
	batch_name = args[0]
	batch_folder = args[1]

	N_bs = 500
	dryrun = False
	verbose = True
	parallel = True
	numprocs = 8

	obs_data = DataReader(batch_name, batch_folder, load_file=None, exclude_fobs=["topct"], dryrun=dryrun, verbose=verbose)
	# fpath = obs_data.write_single_file()
	# obs_data = DataReader(directory_tree)
	# obs_data.retrieve_observable_data(exclude_fobs=["topct"])


	# TEMP TEST
	if (batch_name == "beta60") or (batch_name == "beta6_0"):
		NCfgs = 1000
	else:
		NCfgs = 500

	# # fpath = "../DataGiovanni/24_6.00.npy"
	# # fpath = "../DataGiovanni/28_6.10.npy"
	# fpath = "../DataGiovanni/32_6.20.npy"
	# NCfgs = 500
	# batch_folder = "DataGiovanni"
	# batch_name = "beta62"
	# # fpath = "../data5/beta60/24_6.00.npy"
	# obs_data = DataReader(batch_name, batch_folder, load_file=fpath, NCfgs=NCfgs, exclude_fobs=["topct"], dryrun=dryrun, verbose=verbose)

	# Analysis timers
	pre_time = time.clock()
	observable_strings = []

	# Analyses plaquette data if present in arguments
	if 'plaq' in args:
		plaq_analysis = AnalysePlaquette(obs_data("plaq"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
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
		plaq_analysis.plot_histogram(0)
		plaq_analysis.plot_histogram(-1)
		plaq_analysis.plot_integrated_correlation_time()
		plaq_analysis.plot_integrated_correlation_time()
		plaq_analysis.save_post_analysis_data()

	if 'topc' in args:
		topc_analysis = AnalyseTopologicalCharge(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
		if topc_analysis.beta == 6.0:
			topc_analysis.y_limits = [-9, 9]
		elif topc_analysis.beta == 6.1:
			topc_analysis.y_limits = [-12, 12]
		elif topc_analysis.beta == 6.2:
			topc_analysis.y_limits = [-12, 12]
		else:
			topc_analysis.y_limits = [None, None]

		topc_analysis.boot(N_bs)
		topc_analysis.jackknife()
		topc_analysis.plot_original()
		topc_analysis.plot_boot()
		topc_analysis.plot_jackknife()
		topc_analysis.autocorrelation()
		topc_analysis.plot_autocorrelation(0)
		topc_analysis.plot_autocorrelation(-1)
		topc_analysis.plot_mc_history(0)
		topc_analysis.plot_mc_history(-1)
		topc_analysis.plot_boot()
		topc_analysis.plot_original()
		topc_analysis.plot_jackknife()
		topc_analysis.plot_histogram(0)
		topc_analysis.plot_histogram(-1)
		topc_analysis.plot_integrated_correlation_time()
		topc_analysis.plot_integrated_correlation_time()
		topc_analysis.save_post_analysis_data()

	if 'topcq4' in args:
		topcq4_analysis = AnalyseQQuartic(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
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
		topcq4_analysis.plot_histogram(0)
		topcq4_analysis.plot_histogram(-1)
		topcq4_analysis.plot_integrated_correlation_time()
		topcq4_analysis.plot_integrated_correlation_time()

	if 'qtqzero' in args:
		qzero_flow_time = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99] # Percents of data where we do qtq0
		qtqzero_analysis = AnalyseQtQZero(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
		for qzero_flow_time_index in qzero_flow_time:
			# qtqzero_analysis.y_limits = [0, 2]
			qtqzero_analysis.setQ0(qzero_flow_time_index, y_label=r"$\langle Q_{t}Q_{t_0} \rangle^{1/4} [GeV]$")

			qtqzero_analysis.boot(N_bs)
			qtqzero_analysis.jackknife()
			qtqzero_analysis.plot_original()
			qtqzero_analysis.plot_jackknife()
			qtqzero_analysis.plot_boot()
			# qtqzero_analysis.boot(N_bs, F=qtqzero_analysis.chi, F_error=qtqzero_analysis.chi_std)
			# qtqzero_analysis.jackknife(F=qtqzero_analysis.chi, F_error=qtqzero_analysis.chi_std)
			# qtqzero_analysis.plot_original(correction_function=lambda x: np.power(x, 0.25))
			# qtqzero_analysis.plot_jackknife(correction_function=lambda x: np.power(x, 0.25))
			# qtqzero_analysis.plot_boot(correction_function=lambda x: np.power(x, 0.25))

	if 'topsus' in args:
		topsus_analysis = AnalyseTopologicalSusceptibility(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
		topsus_analysis.boot(N_bs, F=topsus_analysis.chi, F_error=topsus_analysis.chi_std, store_raw_bs_values=True)
		topsus_analysis.jackknife(F=topsus_analysis.chi, F_error=topsus_analysis.chi_std, store_raw_jk_values=True)
		# topsus_analysis.y_limits  = [0.05,0.5]
		topsus_analysis.plot_original()
		topsus_analysis.plot_boot()
		topsus_analysis.plot_jackknife()

		# print "exitsz @ 1096"
		# exit(1)

		topsus_analysis.autocorrelation()
		topsus_analysis.plot_autocorrelation(0)
		topsus_analysis.plot_autocorrelation(-1)
		topsus_analysis.plot_mc_history(0)
		topsus_analysis.plot_mc_history(-1)
		topsus_analysis.plot_original()
		topsus_analysis.plot_boot()
		topsus_analysis.plot_jackknife()
		topsus_analysis.plot_histogram(0)
		topsus_analysis.plot_histogram(-1)
		topsus_analysis.plot_integrated_correlation_time()
		topsus_analysis.plot_integrated_correlation_time()
		topsus_analysis.save_post_analysis_data()

	if 'energy' in args:
		r0 = 0.5
		energy_analysis = AnalyseEnergy(obs_data("energy"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
		energy_analysis.boot(N_bs, store_raw_bs_values=True)
		energy_analysis.jackknife(store_raw_jk_values=True)

		# In case we are at data5, and the x values has already been scaled
		if batch_folder == "data5":
			x_values = energy_analysis.x / r0**2 * energy_analysis.a**2
		else:
			energy_analysis.correction_function_factor = energy_analysis.flow_epsilon ** 2
			x_values = energy_analysis.flow_epsilon * energy_analysis.x / r0**2 * energy_analysis.a**2

		energy_analysis.plot_original(x=x_values, correction_function=energy_analysis.correction_function)
		energy_analysis.plot_boot(x=x_values, correction_function=energy_analysis.correction_function)
		energy_analysis.plot_jackknife(x=x_values, correction_function=energy_analysis.correction_function)
		energy_analysis.autocorrelation()
		energy_analysis.plot_autocorrelation(0)
		energy_analysis.plot_autocorrelation(-1)
		energy_analysis.plot_mc_history(0)
		energy_analysis.plot_mc_history(-1)
		energy_analysis.plot_original(x=x_values, correction_function=energy_analysis.correction_function)
		energy_analysis.plot_boot(x=x_values, correction_function=energy_analysis.correction_function)
		energy_analysis.plot_jackknife(x=x_values, correction_function=energy_analysis.correction_function)
		energy_analysis.plot_histogram(0)
		energy_analysis.plot_histogram(-1, x_label=r"$E$[GeV]")
		energy_analysis.plot_integrated_correlation_time()
		energy_analysis.plot_integrated_correlation_time()
		energy_analysis.save_post_analysis_data()

	if 'topct' in args:
		topct_analysis = AnalyseTopologicalChargeInEuclideanTime(obs_data("topct"),
			batch_name, dryrun=dryrun, parallel=parallel, numprocs=numprocs, 
			verbose = verbose)
		# topct_analysis.boot(N_bs, store_raw_bs_values = True)
		topct_analysis.jackknife(store_raw_jk_values=True)

		topct_analysis.plot_jackknife(x=range(topct_analysis.N_observables_per_config))
		# topct_analysis.plot_original(x = range(topct_analysis.N_observables_per_config))
		# topct_analysis.plot_boot(x = range(topct_analysis.N_observables_per_config))

	post_time = time.clock()
	print "="*100
	print "Analysis of batch %s observables %s in %.2f seconds" % (batch_name, ", ".join([i.lower() for i in args[2:]]), (post_time-pre_time))
	print "="*100

if __name__ == '__main__':
	if not sys.argv[1:]:
		# args = [['beta6_0','data2','plaq','topc','energy','topsus','qtqzero','topcq4','topcqq'],
		# 		['beta6_1','data2','plaq','topc','energy','topsus','qtqzero','topcq4','topcqq'],
		# 		['beta6_2','data2','plaq','topc','energy','topsus','qtqzero','topcq4','topcqq']]

		# args = [['beta6_1','data4','qtqzero']]

		args = [['beta60','data5','topc','plaq','topsus','energy','qtqzero','topcq4'],
				['beta61','data5','topc','plaq','topsus','energy','qtqzero','topcq4']]

		# args = [['beta60', 'data5', 'qtqzero'],
		# 		['beta61', 'data5', 'qtqzero']]

		# args = [['beta6_2','data4','topsus']]
		# args = [['beta61','data5','topc']]
		# args = [['beta61','data5','topc','plaq','topsus','energy','qtqzero','topcq4']]

		args = [['beta60','data5','topc','plaq','topsus','energy','qtqzero']]

		# args = [['test_run_new_counting','output','topc','plaq','energy','topsus']]

		for a in args:
			main(a)
	else:
		main(sys.argv[1:])

	# plt.show()