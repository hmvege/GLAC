import numpy as np, matplotlib.pyplot as plt, sys, os, time

__all__ = ["Autocorrelation","PropagatedAutocorrelation"]

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

			time_used = t2-t1
			args[0].time_used = time_used
			
			print "Autocorrelation: time used with function %s: %.10f secs/ %.10f minutes" % (func.__name__, time_used, time_used/60.)
		
		return val

	return wrapper

class _AutocorrelationCore(object):
	def __init__(self, data, data_std = False, use_numpy = False, data_statistic = lambda x : x, time_autocorrelation = False):
		"""
		Args:
			data 					 (numpy array): dataset to get autocorrelation for
			[optional] use_numpy			(bool): uses numpy to perform the autocorrelation
			[optional] data_statistics	(function):	statistics that should be performed on the data before one proceed with the autocorrelation
			[optional] time_autocorrealtion	(bool): times the autocorrelation function
		Returns:
			Object containing the autocorrelation values
		"""
		# Timer variables
		self.time_autocorrelation = time_autocorrelation
		self.time_used = 0.0

		# Autocorrelation variables
		self.N = len(data)
		self.data = data_statistic(data)
		if data_std != False:
			self.data_std = data_std
		else:
			self.data_std = False
		self.C0 = np.var(data)
		self.R = np.zeros(self.N/2)
		self.R_error = np.zeros(self.N/2)
		self.tau_int = 0
		self.tau_int_error = 0

		# Gets the autocorrelations
		if use_numpy:
			self._numpy_autocorrelation(self.data,self.data)
		else:
			self._autocorrelation(self.data,self.data)

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
		for t in xrange(0,self.N/2):
			for i in xrange(0, self.N - t):
				self.R[t] += (data_x[i] - avg_data_x)*(data_y[i+t] - avg_data_y)
				self.G[t]
			self.R[t] /= (self.N - t)

		self.R /= self.C0

	@timing_function
	def _numpy_autocorrelation(self, data_x, data_y):
		"""
		Numpy method for finding autocorrelation in a dataset. t is tte lag.
		Args:
			Data, (numpy array): dataset to find the autocorrelation in
		Returns:
			C(t)  (numpy array): normalized autocorrelation times 
		"""
		for t in range(0, self.N/2):
			self.R[t] = np.corrcoef(np.array([data_x[0:self.N-t],data_y[t:self.N]]))[0,1]

class Autocorrelation(_AutocorrelationCore):
	"""
	Class for performing an autocorrelation analysis based on Luscher
	"""
	# def __init__(self, data, use_numpy = False, data_statistic = lambda x : x, time_autocorrelation = False):
	def __init__(self,*args,**kwargs):
		"""
		Args:
			data 					 (numpy array): dataset to get autocorrelation for
			[optional] use_numpy			(bool): uses numpy to perform the autocorrelation
			[optional] data_statistics	(function):	statistics that should be performed on the data before one proceed with the autocorrelation
			[optional] time_autocorrealtion	(bool): times the autocorrelation function
		Returns:
			Object containing the autocorrelation values
		"""
		# Calls parent
		super(Autocorrelation,self).__init__(*args,**kwargs)

		# Lambda cutoff
		self.LAMBDA = 100 # As in paper

		# Gets the autocorrelation errors
		map(self._autocorrelation_error,range(self.N/2))

	def integrated_autocorrelation_time(self):
		"""
		Finds the integrated autocorrelation time, and returns it in order to correct the standard deviation
		Returns:
			2*tau_int (float)
		"""
		if self.R[-1] == 0 or self.R[0] == 0:
			print "Error: autocorrelation has not been performed yet!"

		self._get_optimal_w()

		# Sums the integrated autocorrelation time, eq E.12 Luscher(2004)
		self.tau_int = 0.5 + np.sum(self.R[1:self.W])
		return self.tau_int

	def _get_optimal_w(self):
		"""
		Equation E.13 in Luscher(2004)
		"""
		for t in xrange(1,self.N/2):
			if np.abs(self.R[t]) <= np.sqrt(self.R_error[t]):
				self.W = t
				# print self.W, t, self.R[t], np.sqrt(self.R_error[t])
				break
		else:
			self.W = float("NaN")

	def integrated_autocorrelation_time_error(self):
		"""
		Equation E.14 in Luscher(2004)
		"""
		self.tau_int_error = (4*self.W + 2)/float(self.N) * self.tau_int**2
		return self.tau_int_error

	def _autocorrelation_error(self,t):
		"""
		Function for calculating the autocorrelation error.
		Equation E.11 in Luscher
		Args:
			R 		(numpy array): Array of autocorrelations
		Returns:
			R_error (numpy array): Array of error related to the autocorrelation
		"""
		for k in xrange(1,self.LAMBDA + t):
			if k+t >= self.N/2:
				self.R_error[t] += 0.0
			else:
				self.R_error[t] += (self.R[k+t] + self.R[np.abs(k-t)] - 2*self.R[k]*self.R[t])**2 # abs since gamma(t) = gamma(-t)
			# if k-t < 0 or k+t > self.N/2:
			# 	# print "Error: out of bounds at t = %d, k = %d" % (t,k)
			# 	# self.R_error[t] = float("NaN")
			# 	break
		self.R_error[t] /= float(self.N)
		return self.R_error[t]

	def plot_autocorrelation(self, title, filename, lims = 1,dryrun=False):
		"""
		Plots the autocorrelation.
		"""
		fig = plt.figure(dpi=200)
		ax = fig.add_subplot(111)
		# ax.plot(range(self.N/2),self.R,color="r",label="Autocorrelation")
		ax.errorbar(range(self.N/2),self.R,yerr=self.R_error,color="0",ecolor="r",label="Autocorrelation")
		ax.set_ylim(-lims,lims)
		ax.set_xlim(0,self.N/2)
		ax.set_xlabel(r"Lag $t$")
		ax.set_ylabel(r"$R = \frac{C_t}{C_0}$")
		ax.set_title(title,fontsize=16)
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.2))
		ax.grid(True)
		ax.legend()
		if dryrun:
			fig.savefig("tests/autocorrelation_%s.png" % filename)

class PropagatedAutocorrelation(_AutocorrelationCore):
	"""
	Class for performing a general autocorrelation analysis according to 2006 paper by Wolff.

	Assumptions throughout the program:
	- only have 1 alpha, that is only one observable. This simplifies quite alot.
	"""
	# def __init__(self, data, function, function_derivative, use_numpy = False, time_autocorrelation = False):
	def __init__(self, *args, **kwargs):
		"""
		Args:
			data 			 		 (numpy array): dataset to get autocorrelation for
			function 					(function): function to return data as
			function_derivative 		(function):	function derivative to propagate error through
			[optional] use_numpy			(bool): uses numpy to perform the autocorrelation			
			[optional] time_autocorrealtion	(bool): times the autocorrelation function
		Returns:
			Object containing the autocorrelation values
		"""
		# Retrieves relevant functions for later
		self.function = kwargs["function"]
		self.function_derivative = kwargs["function_derivative"]
		del kwargs["function"]
		del kwargs["function_derivative"]

		# Calls parent
		super(PropagatedAutocorrelation,self).__init__(*args,**kwargs)

		# Gets the autocorrelation errors
		self._ac_error();
		# np.asarray(map(self._autocorrelation_error,range(self.N/2)))

	def _ac_error(self):
		# Eq. 6, 7
		avg_data = np.mean(self.data)

		# Eq. 14
		derfun_avg = self.function_derivative(avg_data)

		# Eq. 33
		for t in xrange(self.N/2):
			self.R[t] *= derfun_avg**2

		# Eq. 35, array with different integration cutoffs
		CfW = np.zeros(self.N/2 - 1)
		for W in xrange(1,self.N/2):
			CfW[W-1] = self.R[0] + 2.0*np.sum(self.R[1:W+1])

		# Eq. 49, bias correction
		CfW = CfW*((2*np.arange(len(CfW)) + 1.0)/float(self.N) + 1.0)
		
		# print CfW
		# sys.exit("Exits at 265")

		# Eq. 34
		sigma0 = CfW[0]

		# Eq. 41
		self.tau_int = CfW / (2*sigma0)
		
		# print tau_int
		# sys.exit("Exits at 276")

		for S in np.linspace(1.0,2.0,100):

			# Automatic windowing
			tau = []
			for it, itauint, in enumerate(tau_int):
				if itauint <= 0.5:
					tau.append(0.00000001)
				else:
					# Eq. 51
					tau.append(S/np.log((2*itauint + 1) / (2*itauint - 1)))
			tau = np.asarray(tau)

			# print tau
			# sys.exit("Exits at 288")


			# Eq. 52, getting optimal W
			gW = lambda W,tau : np.exp(-W/tau) - tau/np.sqrt(W*float(self.N))
			for iW, itau in enumerate(tau):
				if iW == 0: continue
				if gW(iW,itau) < 0.0:
					WOptimal = iW
					break
			else:
				WOptimal = float('NaN')

			print WOptimal, S

		# Eq. 42, getting tauint variance
		self.tau_int_error = np.asarray([np.sqrt(4/float(N)*(float(iW) + 0.5 - itau)*itau**2) for iW, itau in tau_int])

		self.tau_int_optimal = self.tau_int[WOptimal]
		self.tau_int_optimal_error = self.tau_int_error[WOptimal]

	def _autocorrelation_error(self,t):
		pass

	def integrated_autocorrelation_time_error(self):
		pass

	def integrated_autocorrelation_time(self):
		pass

def testRegularAC(data,N_bins,store_plots,time_ac_functions):
	print "="*20, "RUNNING DEFAULT TEST", "="*20
	
	# Autocorrelation
	ac = Autocorrelation(data,time_autocorrelation = time_ac_functions)
	ac.plot_autocorrelation(r"Autocorrelation for Plaquette $\beta = 6.2, \tau=10.0$", "beta6_2",dryrun=(not store_plots))
	ac.integrated_autocorrelation_time()
	ac.integrated_autocorrelation_time_error()

	# Autocorrelation with numpy
	ac_numpy = Autocorrelation(data,use_numpy=True, time_autocorrelation = time_ac_functions)
	ac_numpy.plot_autocorrelation(r"Autocorrelation for Plaquette $\beta = 6.2, \tau=10.0$ using Numpy", "beta6_2",dryrun=(not store_plots))
	ac_autocorr = ac_numpy.integrated_autocorrelation_time()
	ac_autocorr_err = ac_numpy.integrated_autocorrelation_time_error()
	print """
Plaquette
Average:                        {0:<.8f}
Std:                            {1:<.8f}
Std with ac-time correction:    {2:<.8f}
sqrt(2*tau_int):                {3:<.8f}
Integrated ac-time:             {4:<.8f}
Integrated ac-time error:       {5:<.8f}""".format(
	np.average(data),
	np.std(data),
	np.std(data)*np.sqrt(2*ac_autocorr),
	np.sqrt(2*ac_autocorr),
	ac_autocorr,
	ac_autocorr_err)

	# Differences in value
	print """
Time used by default method:    {0:<.8f}
Time used by numpy:             {1:<.8f}
Improvement(default/numpy):     {2:<.3f}""".format(ac.time_used, ac_numpy.time_used, ac.time_used/ac_numpy.time_used)
	fig = plt.figure(dpi=200)
	ax = fig.add_subplot(111)
	ax.semilogy(np.abs(ac.R - ac_numpy.R))
	ax.set_title("Relative difference between numpy method and standard method", fontsize=14)
	ax.grid(True)
	if store_plots:
		fig.savefig("tests/relative_differences_in_ac_methods.png")

def testFullAC(data,N_bins,store_plots,time_ac_functions):
	print "="*20, "RUNNING FULL AC TEST", "="*20
	
	def chi_beta6_2(Q_squared):
		const = 0.0763234462734
		return const*Q_squared**(0.25)

	def chi_beta6_2_error(Q_squared):
		const = 0.0763234462734
		return 0.25*const / Q_squared**(0.75)

	ac = PropagatedAutocorrelation(data,function=chi_beta6_2,function_derivative=chi_beta6_2_error,use_numpy=True, time_autocorrelation = time_ac_functions)

def main():
	# Data to load and analyse
	data_plaq = np.loadtxt("tests/plaq_beta6_2_t10.dat")
	data_topc = (np.loadtxt("tests/topc_beta6_2_t10.dat"))**2

	# Histogram bins
	N_bins = 20
	store_plots = True
	time_ac_functions = True

	# testRegularAC(data_plaq,N_bins,store_plots,time_ac_functions)
	testFullAC(data_topc,N_bins,store_plots,time_ac_functions)
	plt.show()

if __name__ == '__main__':
	main()