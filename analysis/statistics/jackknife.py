import numpy as np, matplotlib.pyplot as plt, sys, os

__all__ = ["Jackknife"]

class Jackknife:
	"""
	Class for performing a statistical jack knife.
	"""
	def __init__(self, data, F = lambda x: x, jk_statistics = np.average, non_jk_statistics = lambda x : x):
		"""
		Args:
			data 					(numpy array): 	dataset to give
		Returns:
			Object containing jack-knifed values
		"""
		# Sets some global class variables
		self.N = len(data)
		self.shape = self.N

		# Performs jackknife and sets variables
		self.jk_data_raw = np.zeros((self.N,self.N-1)) # Jack knifed data
		for i in xrange(self.N):
			self.jk_data_raw[i] = np.concatenate([data[:i],data[i+1:]])

		self.jk_data = F(jk_statistics(self.jk_data_raw,axis=1))
		# Estimates variance according to new MHJ book
		self.jk_var = np.var(self.jk_data)*(self.N - 1)
		self.jk_std = np.sqrt(self.jk_var)
		self.jk_avg = np.average(self.jk_data)
		
		# Gets and sets non-bootstrapped values
		data = F(non_jk_statistics(data))
		self.avg_original = np.average(data)
		self.var_original = np.var(data)
		self.std_original = np.std(data)

		# Returns the unbiased estimator/average
		self.jk_avg_unbiased = self.avg_original - (self.N - 1) * (self.jk_avg - self.avg_original)

	def __call__(self):
		"""
		When called, returns the bootstrapped samples.
		"""
		return self.jk_data

	def __len__(self):
		"""
		Length given as number of bootstraps.
		"""
		return self.shape[0]

	def __str__(self):
		"""
		When object is printed, prints information about the jackknife performed.
		"""
		msg = """
JACKKNIFE:
%s
Non-jackknife: %10.10f %10.10E %10.10E
Jackknife:     %10.10f %10.10E %10.10E %10.10E(unbiased average)
		""" % ("="*61,self.avg_original, self.var_original, self.std_original, self.jk_avg, self.jk_var, self.jk_std, self.jk_avg_unbiased)
		return msg

def main():
	# Data to load and analyse
	data = np.loadtxt("tests/plaq.dat",skiprows=8)
	
	# Histogram bins
	N_bins = 20
	
	# Jackknifing
	jk = Jackknife(data)
	jk_data = jk()

	print jk

	plt.show()

if __name__ == '__main__':
	main()
