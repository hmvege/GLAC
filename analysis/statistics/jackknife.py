import numpy as np, matplotlib.pyplot as plt, sys, os

__all__ = ["Jackknife"]

class Jackknife:
	"""
	Class for performing a statistical jack knife.
	"""
	def __init__(self, data, F = lambda x: x, data_statistics = np.average): # , statistics = np.mean, F = lambda x : x, statistics = lambda x: x
		"""
		Args:
			data 					(numpy array): 	dataset to give
		Returns:
			Object containing jack-knifed values
		"""
		# Sets some global class variables
		self.N = len(data)
		self.shape = self.N


		# data = F(data_statistics(data))

		# Gets and sets non-bootstrapped values
		self.avg_original = data_statistics(data)
		self.var_original = np.var(data)
		self.std_original = np.std(data)

		# Performs jackknife and sets variables
		self.jk_data = np.zeros(self.N) # Jack knifed data
		for i in xrange(self.N):
			self.jk_data[i] = np.average(np.concatenate([data[:i-1],data[i:]]))

		# F(data_statistics())
		# self.jk_var = self.__calculate_variance(self.jk_data, self.avg_original)
		self.jk_var = np.var(self.jk_data)
		self.jk_std = np.sqrt(self.jk_var)
		self.jk_avg_biased = data_statistics(self.jk_data)

		# Returns the unbiased estimator/average
		self.jk_avg = self.avg_original - (self.N - 1) * (data_statistics(self.jk_data) - self.avg_original)

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

	def get_average(self):
		"""
		Returns the biased average of the jack knife method.
		"""
		return self.jk_avg

	def get_average_unbiased(self):
		"""
		Returns the unbiased average of the jack knife method.
		"""
		return self.jk_avg_unbiased

	def __calculate_variance(self,jackknifed_dataset, avg_orginal):
		"""
		Returns the variance of the jack knife method.
		"""
		var = 0
		for i in xrange(self.N):
			var += (jackknifed_dataset[i] - avg_orginal)**2
		var = (self.N - 1) / float(self.N) * var
		return var

	def __str__(self):
		"""
		When object is printed, prints information about the jackknife performed.
		"""
		msg = """
JACKKNIFE:
%s
Non-jackknife: %10.10f %10.10E %10.10E
Jackknife:     %10.10f %10.10E %10.10E
		""" % ("="*61,self.avg_original, self.var_original, self.std_original, self.jk_avg, self.jk_var, self.jk_std)
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
