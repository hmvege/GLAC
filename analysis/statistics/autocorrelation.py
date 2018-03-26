import numpy as np, matplotlib.pyplot as plt, sys, os, time

__all__ = ["Autocorrelation", "PropagatedAutocorrelation"]

"""
Books:
Quantum Chromo Dynamics on the Lattice, Gattringer
Papers:
Schwarz-preconditioned HMC algorithm for two-flavour lattice QCD, M. Luscher 2004
Monte Carlo errors with less errors, U. Wolff 2006
"""

def timing_function(func):
	"""
	Time function.
	"""
	def wrapper(*args):
		if args[0].time_autocorrelation:
			t1 = time.clock()
		
		val = func(*args)

		if args[0].time_autocorrelation:
			t2 = time.clock()

			time_used = t2 - t1
			args[0].time_used = time_used
			
			print ("Autocorrelation: time used with function %s: %.10f secs/"
				" %.10f minutes" % (func.__name__, time_used, time_used/60.))
		
		return val

	return wrapper

class _AutocorrelationCore(object):
	def __init__(self, data, function_derivative=lambda x: x, function_parameters=None, method="correlate", time_autocorrelation=False):
		"""
		Args:
			data 					 (numpy array): dataset to get autocorrelation for
			function_derivative			(function): the derivative of function to propagate data through
			function_parameters		  (dictionary): dictionary of function derivative parameters
			[optional] method				 (str): method of performing autocorrelation: "corroeff", "correlate", "manual"
			[optional] time_autocorrealtion	(bool): times the autocorrelation function
		Returns:
			Object containing the autocorrelation values
		"""
		# Retrieves relevant functions for later
		self.function_derivative = function_derivative
		self.function_parameters = function_parameters

		# Timer variables
		self.time_autocorrelation = time_autocorrelation
		self.time_used = 0.0

		# Autocorrelation variables
		self.data = data
		self.N = len(self.data)
		self.C0 = np.var(self.data)
		self.R = np.zeros(self.N/2)
		self.G = np.zeros(self.N/2)
		self.R_error = np.zeros(self.N/2)
		self.tau_int = 0
		self.tau_int_error = 0

		# Gets the autocorrelations
		if method == "corrcoef":
			self._numpy_autocorrelation(self.data, self.data)
		elif method == "correlate":
			self._numpy2_autocorrelation(self.data, self.data)
		elif method == "manual":
			self._autocorrelation(self.data, self.data)
		else:
			raise KeyError("Method of autocorrelation not recognized among: corrcoef, correlate, manual")

	def __call__(self):
		"""
		Returns:
			(numpy double array) the auto-correlation
		"""
		return self.R

	def __len__(self):
		"""
		Returns:
			(int) the length of the auto-correlation results.
		"""
		return len(self.R)

	@timing_function
	def _autocorrelation(self, data_x, data_y):
		"""
		Gets the autocorrelation from a dataset.
		Args:
			Data, (numpy array): dataset to find the autocorrelation in
		Returns:
			C(t)  (numpy array): normalized autocorrelation times 
		"""
		avg_data_x = np.average(data_x)
		avg_data_y = np.average(data_y)
		for t in xrange(0, self.N/2):
			for i in xrange(0, self.N - t):
				self.R[t] += (data_x[i] - avg_data_x)*(data_y[i+t] - avg_data_y)

			self.R[t] /= (self.N - t)

		self.G = self.R * 2
		self.R /= self.C0

	@timing_function
	def _numpy_autocorrelation(self, data_x, data_y):
		"""
		Numpy method for finding autocorrelation in a dataset. t is the lag.
		Args:
			Data, (numpy array): dataset to find the autocorrelation in
		Returns:
			C(t)  (numpy array): normalized autocorrelation times 
		"""
		for t in range(0, self.N/2):
			self.R[t] = np.corrcoef(np.array([data_x[0:self.N - t], data_y[t:self.N]]))[0,1]

		self.G = self.R * self.C0 * 2

	@timing_function
	def _numpy2_autocorrelation(self, x, y):
		"""
		http://stackoverflow.com/q/14297012/190597
		http://en.wikipedia.org/wiki/Autocorrelation#Estimation
		"""
		n = len(x)
		# variance = x.var()
		x = x-x.mean()
		y = y-y.mean()
		self.G = np.correlate(x, y, mode="full")[-n:]
		self.G /= np.arange(n, 0, -1)
		self.G = self.G[:self.N/2]
		self.R = self.G/self.G[0]

	def integrated_autocorrelation_time(self, plot_cutoff=False):
		raise NotImplemented("integrated_autocorrelation_time() not implemented for base class")

	def integrated_autocorrelation_time_error(self):
		raise NotImplemented("integrated_autocorrelation_time_error() not implemented for base class")

	def plot_autocorrelation(self, title, filename, lims=1, dryrun=False):
		"""
		Plots the autocorrelation.
		"""
		x = range(len(self.R))
		y = self.R
		y_std = self.R_error

		fig = plt.figure(dpi=200)
		ax = fig.add_subplot(111)
		ax.plot(x, y, color="0", label="Autocorrelation")
		ax.fill_between(x, y - y_std, y + y_std, alpha=0.5, edgecolor='', facecolor='r')
		# ax.errorbar(x,y,yerr=y_std,color="0",ecolor="r",label="Autocorrelation")
		ax.set_ylim(-lims, lims)
		ax.set_xlim(0, self.N/2)
		ax.set_xlabel(r"Lag $t$")
		ax.set_ylabel(r"$R = \frac{C_t}{C_0}$")
		ax.set_title(title, fontsize=16)
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.2))
		ax.grid(True)
		ax.legend()
		if dryrun:
			fig.savefig("tests/autocorrelation_%s.png" % filename)

class Autocorrelation(_AutocorrelationCore):
	"""
	Class for performing an autocorrelation analysis based on Luscher
	"""
	# def __init__(self, data, use_numpy = False, data_statistic = lambda x : x, time_autocorrelation = False):
	def __init__(self, *args, **kwargs):
		# Calls parent
		super(Autocorrelation, self).__init__(*args, **kwargs)

		# Lambda cutoff
		self.LAMBDA = 100 # As in paper

		# Gets the autocorrelation errors
		map(self._autocorrelation_error, range(self.N/2))

	def _autocorrelation_error(self, t):
		"""
		Function for calculating the autocorrelation error.
		Equation E.11 in Luscher
		Args:
			R 		(numpy array): Array of autocorrelations
		Returns:
			R_error (numpy array): Array of error related to the autocorrelation
		"""
		for k in xrange(1, self.LAMBDA + t):
			if k+t >= self.N/2:
				break
			else:
				self.R_error[t] += (self.R[k+t] + self.R[np.abs(k-t)] - 2*self.R[k]*self.R[t])**2 # abs since gamma(t) = gamma(-t)

		self.R_error[t] = np.sqrt(self.R_error[t] / float(self.N))
		return self.R_error[t]

	def _get_optimal_w(self):
		"""
		Equation E.13 in Luscher(2004)
		"""
		for t in xrange(1,self.N/2):
			if np.abs(self.R[t]) <= self.R_error[t]:
				self.W = t
				break
		else:
			self.W = float("NaN")

	def integrated_autocorrelation_time(self, plot_cutoff=False):
		"""
		Finds the integrated autocorrelation time, and returns it in order to correct the standard deviation
		Returns:
			2*tau_int (float)
		"""
		if self.R[-1] == 0 or self.R[0] == 0:
			print "Error: autocorrelation has not been performed yet!"

		self.tau_int = np.array([0.5 + np.sum(self.R[1:iW]) for iW in xrange(1, self.N/2)])

		self._get_optimal_w()

		# Sums the integrated autocorrelation time, eq E.12 Luscher(2004)
		# self.tau_int = 0.5 + np.sum(self.R[1:self.W])

		self.tau_int_optimal = self.tau_int[self.W]

		# Plots cutoff if prompted
		if plot_cutoff:
			tau = np.zeros(len(self.R)-1)
			for iW in xrange(1, len(self.R)):
				tau[iW] = (0.5 + np.sum(self.R[1:iW]))

			plt.figure()
			plt.plot(range(1,len(self.R)),tau)
			plt.title(r"Cutoff $W = %d$" % (self.W))
			plt.xlabel(r"$W$")
			plt.ylabel(r"$\tau_{int}(W)$")
			plt.axvline(self.W)
			plt.show()

		return self.tau_int_optimal

	def integrated_autocorrelation_time_error(self):
		"""
		Equation E.14 in Luscher(2004)
		"""
		self.tau_int_error = np.sqrt((4*self.W + 2)/float(self.N) * self.tau_int**2)
		# self.tau_int_error = np.sqrt(4/float(self.N) * (self.W + 0.5 - self.tau_int) * self.tau_int**2)
		self.tau_int_optimal_error = self.tau_int_error[self.W]
		return self.tau_int_optimal_error

class PropagatedAutocorrelation(_AutocorrelationCore):
	"""
	Class for performing a general autocorrelation analysis according to 2006 paper by Wolff.

	Assumptions throughout the program:
	- only have 1 alpha, that is only one observable. This simplifies quite alot.
	"""
	def __init__(self, *args, **kwargs):
		# Calls parent
		super(PropagatedAutocorrelation, self).__init__(*args, **kwargs)

		# Gets the autocorrelation errors
		self._autocorrelation_error();

	def _autocorrelation_error(self, SParam=1.11):
		# Eq. 6, 7
		avg_data = np.mean(self.data)

		# Eq. 14
		derfun_avg = self.function_derivative(avg_data, **self.function_parameters)

		# Eq. 33
		self.G *= derfun_avg**2

		# Eq. 35, array with different integration cutoffs
		CfW = np.array([self.G[0]] + [self.G[0] + 2.0*np.sum(self.G[1:W+1]) for W in xrange(1, self.N/2)])
		# CfW = np.array([0.5 + np.sum(self.R[1:iW]) for iW in xrange(1,self.N/2)])

		# Eq. 49, bias correction
		CfW = self._correct_bias(CfW)
		
		# Eq. 34
		sigma0 = CfW[0]

		# Eq. 41
		self.tau_int = CfW / (2*sigma0)

		if SParam == False:
			for S in np.linspace(1.0, 2.0, 20):
				self.W = self._automatic_windowing_procedure(S)

				plt.figure()
				plt.plot(range(len(self.tau_int)/2), self.tau_int[:len(self.tau_int)/2])
				plt.title(r"$W = %d$, $S_{param} = %.2f$" % (self.W,S))
				plt.ylim(0,1.25*np.max(self.tau_int[:len(self.tau_int)/4]))
				plt.xlabel(r"$W$")
				plt.ylabel(r"$\tau_{int}(W)$")
				plt.axvline(self.W)
				plt.show()
		else:
			self.W = self._automatic_windowing_procedure(SParam)

		self.tau_int_optimal = self.tau_int[self.W]

	def _correct_bias(self, CfW):
		"""
		Eq. 49, bias correction
		"""
		return CfW*((2*np.arange(len(CfW)) + 1.0)/float(self.N) + 1.0)

	def _gW(self, tau):
		"""
		Eq. 52, getting optimal W
		"""
		for iW, itau in enumerate(tau):
			if iW == 0: continue
			if np.exp(-iW/itau) - itau/np.sqrt(iW*float(self.N)) < 0.0:
				return iW
		else:
			return float('NaN')

	def _automatic_windowing_procedure(self, S):
		"""
		Automatic windowing as described in Wolff paper, section 3.3 
		"""
		tau = []
		for it, itauint, in enumerate(self.tau_int):
			if itauint <= 0.5:
				tau.append(0.00000001)
			else:
				# Eq. 51
				tau.append(S/np.log((2*itauint + 1) / (2*itauint - 1)))
		tau = np.asarray(tau)

		return self._gW(tau)

	def integrated_autocorrelation_time_error(self):
		"""
		Eq. 42, standard deviation of tau_int
		"""
		self.tau_int_error = np.asarray([np.sqrt(4/float(self.N)*(float(iW) + 0.5 - itau)*itau**2) for iW, itau in enumerate(self.tau_int)])
		# self.tau_int_error = np.sqrt((4*self.W + 2)/float(self.N) * self.tau_int**2)
		self.tau_int_optimal_error = self.tau_int_error[self.W]
		return self.tau_int_optimal_error

	def integrated_autocorrelation_time(self):
		return self.tau_int_optimal

def testRegularAC(data, N_bins, store_plots, time_ac_functions):
	"""Function for testing default autocorrelation method."""
	
	def chi_beta6_2_derivative(Q_squared):
		const = 0.0763234462734
		return 0.25*const / Q_squared**(0.75)

	def print_values(observable, method, values, autocorr, autocorr_error):
		value_string = observable
		value_string += "\nMethod:                             {0:<s}".format(method)
		value_string += "\nAverage:                            {0:<.8f}".format(np.average(values))
		value_string += "\nStd:                                {0:<.8f}".format(np.std(values))
		value_string += "\nStd with ac-time correction:        {0:<.8f}".format(np.std(data)*np.sqrt(2*autocorr))
		value_string += "\nsqrt(2*tau_int):                    {0:<.8f}".format(np.sqrt(2*autocorr))
		value_string += "\nIntegrated ac-time:                 {0:<.8f}".format(autocorr)
		value_string += "\nIntegrated ac-time error:           {0:<.8f}".format(autocorr_error)
		print value_string

	print "="*20, "RUNNING DEFAULT TEST", "="*20
	
	# Autocorrelation
	ac = Autocorrelation(data, method="manual", time_autocorrelation=time_ac_functions)
	ac.plot_autocorrelation(r"Autocorrelation for Plaquette $\beta = 6.2, \tau=10.0$", "beta6_2", dryrun=(not store_plots))
	ac_manual = ac.integrated_autocorrelation_time()
	ac_manual_err = ac.integrated_autocorrelation_time_error()

	# Autocorrelation with numpy corrcoef
	ac_numpy1 = Autocorrelation(data, method="corrcoef", time_autocorrelation=time_ac_functions)
	ac_numpy1.plot_autocorrelation(r"Autocorrelation for Plaquette $\beta = 6.2, \tau=10.0$ using np\.corrcoef", "beta6_2", dryrun=(not store_plots))
	ac_autocorr = ac_numpy1.integrated_autocorrelation_time()
	ac_autocorr_err = ac_numpy1.integrated_autocorrelation_time_error()

	ac_numpy2 = Autocorrelation(data, method="correlate", time_autocorrelation=time_ac_functions)
	ac_numpy2.plot_autocorrelation(r"Autocorrelation for Plaquette $\beta = 6.2, \tau=10.0$ using np\.correlate", "beta6_2", dryrun=(not store_plots))
	ac_autocorr2 = ac_numpy2.integrated_autocorrelation_time()
	ac_autocorr_err2 = ac_numpy2.integrated_autocorrelation_time_error()

	print_values("Plaquette", "manual", data, ac_manual, ac_manual_err)
	print_values("Plaquette", "corrcoef", data, ac_autocorr, ac_autocorr_err)
	print_values("Plaquette", "correlate", data, ac_autocorr2, ac_autocorr_err2)

	# Differences in time
	print """
Time used by default method:    	{0:<.8f}
Time used by numpy corrcoef:    	{1:<.8f}
Time used by numpy correlate:   	{2:<.8f}
Improvement(default/corrcoef): 		{3:<.3f}
Improvement(default/correlate): 	{4:<.3f}
Improvement(corrcoef/correlate):	{5:<.3f}""".format(ac.time_used,ac_numpy1.time_used, ac_numpy2.time_used, ac.time_used/ac_numpy1.time_used, ac.time_used/ac_numpy2.time_used,ac_numpy1.time_used/ac_numpy2.time_used)

	# Plotting difference
	fig = plt.figure(dpi=200)
	ax = fig.add_subplot(111)
	ax.semilogy(np.abs(ac.R - ac_numpy1.R))
	ax.set_title("Relative difference between numpy method and standard method", fontsize=14)
	ax.grid(True)
	if store_plots:
		fig.savefig("tests/relative_differences_in_ac_methods.png")

def testFullAC(data,N_bins,store_plots,time_ac_functions):
	"""Function for testing autocorrelation with error propagation."""

	print "="*20, "RUNNING FULL AC TEST", "="*20

	def chi_beta6_2_derivative(Q_squared):
		const = 0.0763234462734
		return 0.25*const / Q_squared**(0.75)

	ac1 = Autocorrelation(data,method="corrcoef",time_autocorrelation=time_ac_functions)
	ac1.plot_autocorrelation(r"Autocorrelation for Topological Suscpetibility $\beta = 6.2$","beta6_2_topc",dryrun=(not store_plots))
	print ac1.integrated_autocorrelation_time()
	print ac1.integrated_autocorrelation_time_error()
	print ac1.W

	# ac = PropagatedAutocorrelation(data,function_derivative=chi_beta6_2_derivative,method="corrcoef",time_autocorrelation=time_ac_functions)
	ac = PropagatedAutocorrelation(data,function_derivative=chi_beta6_2_derivative,method="corrcoef",time_autocorrelation=time_ac_functions)
	ac.plot_autocorrelation(r"Autocorrelation for Topological Suscpetibility $\beta = 6.2$","beta6_2_topc",dryrun=(not store_plots))
	print ac.integrated_autocorrelation_time()
	print ac.integrated_autocorrelation_time_error()
	print ac.W

def main():
	# Data to load and analyse
	data_plaq = np.loadtxt("tests/plaq_beta6_2_t10.dat")
	data_topc = (np.loadtxt("tests/topc_beta6_2_t10.dat"))**2

	# Histogram bins
	N_bins = 20
	store_plots = True
	time_ac_functions = True

	# #### THESE METHODS SHOULD GIVE EQUAL ANSWERS!! --> and they now mostly do
	# ac1 = Autocorrelation(data_plaq,use_numpy=False,time_autocorrelation=time_ac_functions)
	# ac = PropagatedAutocorrelation(data_plaq,use_numpy=False,function=lambda x:x,function_derivative=lambda x:x,time_autocorrelation=time_ac_functions)
	# print "{0:<15s} {1:<15s} {2:<5s}".format("tau_int","tau_int_eror","w")
	# p = lambda t,terr,w : "{0:<15.10f} {1:<15.10f} {2:<5d}".format(t,terr,w)
	# print p(ac1.integrated_autocorrelation_time(),ac1.integrated_autocorrelation_time_error(),ac1.W)	
	# print p(ac.integrated_autocorrelation_time(),ac.integrated_autocorrelation_time_error(),ac.W)

	# testRegularAC(data_plaq,N_bins,store_plots,time_ac_functions)
	testRegularAC(data_plaq,N_bins,store_plots,time_ac_functions)
	testFullAC(data_topc,N_bins,store_plots,time_ac_functions)
	# plt.show()

if __name__ == '__main__':
	main()