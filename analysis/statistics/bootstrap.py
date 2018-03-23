import numpy as np, matplotlib.pyplot as plt, sys, os

__all__ = ["Bootstrap"]

class Bootstrap:
	"""
	Class for creating a bootstrap sample.
	"""
	def __init__(self, data, N_BS, index_lists=[], seed=None, axis=None):
		"""
		Args:
			data 					(numpy array): 	dataset to give
			N_BS					(int): 			number of bootstrap-samples to create
			*bootstrap_statistics	(function):		statistics to run on bootstrap sample, default is numpy.mean()
			*index_lists			(numpy array): 	randomly generated lists that can be provided
			*seed 					(int/float):	seeding the random list generation
		Returns:
			Object containing bootstrapped values
		"""
		N = len(data)
		if seed != None: # Generates a seed if it is provided
			np.random.seed(seed=seed)
			self.seed = seed
		if len(index_lists) == 0: # Allows user to send in a predefined list if needed
			index_lists = np.random.randint(N, size=(N_BS, N))
		self.bs_data_raw = data[index_lists]
		self.bs_data = np.mean(self.bs_data_raw, axis=1)

		# Performing basic bootstrap statistics
		self.bs_avg = np.average(self.bs_data, axis=axis)
		self.bs_var = np.var(self.bs_data, axis=axis)
		self.bs_std = np.sqrt(self.bs_var)

		# Performing basic statistics on original data
		self.data_original = data
		self.avg_original = np.average(self.data_original, axis=axis)
		self.var_original = np.var(self.data_original, axis=axis)
		self.std_original = np.std(self.data_original, axis=axis)

		# Sets some global class variables
		self.shape = self.bs_avg.shape
		self.N_BS = N_BS

	def __add__(self, other):
		"""
		Enables adding of two bootstrapped datasets.
		"""
		if type(other) != Bootstrap:
			raise TypeError("%s should be of type Bootstrap." % other)
 		new_data = np.concatenate((self.bs_avg, other.bs_avg), axis=0)
		return new_data

	def __pow__(self, other):
		return np.power(self.bs_avg, other)

	def __call__(self):
		"""
		When called, returns the bootstrapped samples
		"""
		return self.bs_data

	def __len__(self):
		"""
		Length given as number of bootstraps
		"""
		return self.shape[0]

	def __str__(self):
		"""
		When object is printed, prints information about the bootstrap performed.
		"""
		msg = """
BOOTSTRAP:
%s
Non-bootstrap: %10.10f %10.10E %10.10E
Bootstrap:     %10.10f %10.10E %10.10E
N boostraps:   %d
		""" % ("="*61, self.avg_original, self.var_original, self.std_original, self.bs_avg, self.bs_var, self.bs_std, self.N_BS)
		return msg

def main():
	# Data to load and analyse
	data = np.loadtxt("tests/plaq_beta6_2_t10.dat", skiprows=8)

	# Histogram bins
	N_bins = 20

	# Bootstrapping
	N_bootstraps = int(500)
	bs = Bootstrap(data, N_bootstraps)
	bs_data = bs()

	print bs

if __name__ == '__main__':
	main()