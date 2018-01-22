import numpy as np, matplotlib.pyplot as plt, sys, os, time

__all__ = ["Autocorrelation"]

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

class Autocorrelation:
	"""
	Class for performing an autocorrelation analysis.
	"""
	def __init__(self, data, use_numpy = False, time_autocorrelation = False):
		"""
		Args:
			data 							(numpy array): 	dataset to get autocorrelation for
			[optional] use_numpy			(bool): 		uses numpy to perform the autocorrelation
			[optional] time_autocorrealtion	(bool): 		times the autocorrelation function
		Returns:
			Object containing the autocorrelation values
		"""
		# Timer variables
		self.time_autocorrelation = time_autocorrelation
		self.time_used = 0.0

		# Autocorrelation variables
		self.N = len(data)
		self.data = data
		self.C0 = np.var(data)
		if use_numpy:
			self.R = self._get_numpy_autocorrelation(data)
		else:
			self.R = self._get_autocorrelation(data)

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

	def integrated_autocorrelation_time(self):
		"""
		Finds the integrated autocorrelation time, and returns it in order to correct the standard deviation
		Returns:
			2*tau_int (float)
		"""
		self.tau_int = 0.5 + np.sum(self.R)
		return self.tau_int

	@timing_function
	def _get_autocorrelation(self, data):
		"""
		Gets the autocorrelation from a dataset.
		Args:
			Data, (numpy array): dataset to find the autocorrelation in
		Returns:
			C(t)  (numpy array): normalized autocorrelation times 
		"""
		avg_data = np.average(data)
		C = np.zeros(self.N/2)
		for h in xrange((self.N)/2): # PARALLELIZE THIS FOR LOOP!
			for i in xrange(0, self.N - h):
				C[h] += (data[i] - avg_data)*(data[i+h] - avg_data)
			C[h] /= (self.N - h)
		return C / self.C0

	@timing_function
	def _get_numpy_autocorrelation(self, data):
		"""
		Numpy method for finding autocorrelation in a dataset. h is the lag.
		Args:
			Data, (numpy array): dataset to find the autocorrelation in
		Returns:
			C(t)  (numpy array): normalized autocorrelation times 
		"""
		R = np.zeros(self.N/2)
		for h in range(0, self.N/2): # PARALLELIZE THIS FOR LOOP!
			R[h] = np.corrcoef(np.array([data[0:self.N-h],data[h:self.N]]))[0,1]
		return R

	def plot_autocorrelation(self, title, filename, lims = 1,dryrun=False):
		"""
		Plots the autocorrelation.
		"""
		fig = plt.figure(dpi=200)
		ax = fig.add_subplot(111)
		ax.plot(range(self.N/2),self.R,color="r",label="Autocorrelation")
		ax.set_ylim(-lims,lims)
		ax.set_xlim(0,self.N/2)
		ax.set_xlabel(r"Lag $h$")
		ax.set_ylabel(r"$R = \frac{C_h}{C_0}$")
		ax.set_title(title,fontsize=16)
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.2))
		ax.grid(True)
		ax.legend()
		if dryrun:
			fig.savefig("tests/autocorrelation_%s.png" % filename)

# class IntegratedAutoCorrelationTime:
# 	"""
# 	Finds the integrated autocorrelation time
# 	"""
# 	def __init__(self, ac_objects):
# 		"""
# 		Takes an list/array of autocorrelation objects, should be of size NFlows * NConfigurations/2
# 		"""
# 		None

def main():
	# Data to load and analyse
	data = np.loadtxt("tests/plaq.dat",skiprows=8)
	
	# Histogram bins
	N_bins = 20

	store_plots = True
	time_ac_functions = True
	
	# Autocorrelation
	ac = Autocorrelation(data,time_autocorrelation = time_ac_functions)
	ac.plot_autocorrelation(r"Autocorrelation for $\beta = 6.1$", "beta6_1",dryrun=(not store_plots))

	# Autocorrelation with numpy
	ac_numpy = Autocorrelation(data,use_numpy=True, time_autocorrelation = time_ac_functions)
	ac_numpy.plot_autocorrelation(r"Autocorrelation for $\beta = 6.1$ using Numpy", "beta6_1",dryrun=(not store_plots))
	
	# Differences in value
	print "Time used by default method: {0:<.8f}\nTime used by numpy: {1:<.8f}\nImprovement(default/numpy): {2:<.3f}".format(ac.time_used, ac_numpy.time_used, ac.time_used/ac_numpy.time_used)
	fig = plt.figure(dpi=200)
	ax = fig.add_subplot(111)
	ax.semilogy(np.abs(ac.R - ac_numpy.R))
	ax.set_title("Relative difference between numpy method and standard method", fontsize=14)
	ax.grid(True)
	if store_plots:
		fig.savefig("tests/relative_differences_in_ac_methods.png")

	plt.show()

if __name__ == '__main__':
	main()