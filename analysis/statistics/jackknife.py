import numpy as np, matplotlib.pyplot as plt, sys, os, time

__all__ = ["Jackknife"]

def timing_function(func):
	"""
	Time function.
	"""
	def wrapper(*args):
		if args[0].time_jk:
			t1 = time.clock()
		
		val = func(*args)

		if args[0].time_jk:
			t2 = time.clock()

			time_used = t2-t1
			args[0].time_used = time_used
			
			print "Autocorrelation: time used with function %s: %.10f secs/ %.10f minutes" % (func.__name__, time_used, time_used/60.)
		
		return val

	return wrapper

class Jackknife:
	"""
	Class for performing a statistical jack knife.
	"""
	def __init__(self, data, time_jk=False, axis=None):
		"""
		Method for jack-knifing a dataset.

		Args:
			data: data containes in a numpy array
		"""
		# Timer variable
		self.time_jk = time_jk

		# Sets some global class variables
		if len(data.shape) == 2:
			self.N, self.N_array_points = data.shape
			self.jk_data_raw = np.zeros((self.N, self.N-1, self.N_array_points)) # Jack knifed data
		else:
			self.N = len(data)
			self.jk_data_raw = np.zeros((self.N, self.N-1)) # Jack knifed data

		# Performs jackknife and sets variables
		self._perform_jk(data)

		# Performing statistics on jackknife samples
		self.jk_var = np.var(self.jk_data, axis=axis)*(self.N - 1) # Estimates variance according to new MHJ book
		self.jk_std = np.sqrt(self.jk_var) # Ensures that the errors are run through the function F
		self.jk_avg_biased = np.average(self.jk_data, axis=axis)

		# Gets and sets non-bootstrapped values
		self.avg_original = np.average(data, axis=axis)
		self.var_original = np.var(data, axis=axis) # Ensures that the errors are run through the function F
		self.std_original = np.sqrt(self.var_original)

		# Returns the unbiased estimator/average
		self.jk_avg = self.N*self.avg_original - (self.N - 1) * self.jk_avg_biased
		self.jk_avg_unbiased = self.jk_avg

	@timing_function
	def _perform_jk(self,data):
		"""
		Function for performing the jackknife.
		"""
		for i in xrange(self.N):
			# self.jk_data_raw[i] = np.concatenate([data[:i],data[i+1:]])
			self.jk_data_raw[i][0:i] = data[:i]
			self.jk_data_raw[i][i:] = data[i+1:]
		self.jk_data = np.average(self.jk_data_raw, axis=1)

	def __call__(self):
		"""
		When called, returns the bootstrapped samples.
		"""
		return self.jk_data

	def __len__(self):
		"""
		Length given as number of bootstraps.
		"""
		return self.N

	def __str__(self):
		"""
		When object is printed, prints information about the jackknife performed.
		"""
		msg = """
JACKKNIFE:
%s
Non-jackknife: %10.10f %10.10E %10.10E
Jackknife:     %10.10f %10.10E %10.10E %20.16f (unbiased average)
		""" % ("="*61, self.avg_original, self.var_original, self.std_original, self.jk_avg, self.jk_var, self.jk_std, self.jk_avg_unbiased)
		return msg

def main():
	# Data to load and analyse
	data = np.loadtxt("tests/plaq_beta6_2_t10.dat", skiprows=8)

	# Jackknifing
	jk = Jackknife(data, time_jk=True)
	jk_data = jk()

	print jk

if __name__ == '__main__':
	main()
