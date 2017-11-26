import numpy as np, matplotlib.pyplot as plt, sys, os

class Bootstrap:
	"""
	Class for creating a bootstrap sample.
	"""
	def __init__(self, data, N_BS, bootstrap_statistics = np.mean, F = lambda x : x, non_bs_stats = lambda x: x, index_lists=[], seed=None):
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

		# Gets and sets bootstrapped values 
		self.bs_data_raw = data[index_lists]
		self.bs_data = bootstrap_statistics(self.bs_data_raw,axis=1)
		self.bs_data = F(self.bs_data)
		self.bs_avg = np.mean(self.bs_data)
		self.bs_var = np.var(self.bs_data)
		self.bs_std = np.std(self.bs_data)

		# Gets and sets non-bootstrapped values
		data = F(non_bs_stats(data))
		self.avg_original = np.average(data)
		self.var_original = np.var(data)
		self.std_original = np.std(data)

		# Sets some global class variables
		self.shape = self.bs_avg.shape
		self.N_BS = N_BS

	def __add__(self, other):
		"""
		Enables adding of two bootstrapped datasets.
		"""
		if type(other) != Bootstrap:
			raise TypeError("%s should be of type Bootstrap." % other)
 		new_data = np.concatenate((self.bs_avg,other.bs_avg),axis=0)
		return new_data

	def __pow__(self, other):
		return np.power(self.bs_avg,other)

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
		""" % ("="*61,self.avg_original, self.var_original, self.std_original, self.bs_avg, self.bs_var, self.bs_std, self.N_BS)
		return msg

class Jackknife:
	"""
	Class for performing a statistical jack knife.
	"""
	def __init__(self, data):
		"""
		Args:
			data 					(numpy array): 	dataset to give
		Returns:
			Object containing jack-knifed values
		"""
		# Sets some global class variables
		self.N = len(data)
		self.shape = self.N

		# Gets and sets non-bootstrapped values
		self.avg_original = np.average(data)
		self.var_original = np.var(data)
		self.std_original = np.std(data)

		# Performs jackknife and sets variables
		self.jk_data = np.zeros(self.N) # Jack knifed data
		for i in xrange(self.N):
			self.jk_data[i] = np.average(np.concatenate([data[:i-1],data[i:]]))
		self.jk_avg = np.average(self.jk_data)
		self.jk_avg_unbiased = self.avg_original - (self.N - 1) * (self.jk_avg - self.avg_original)
		self.jk_var = np.var(self.jk_data)
		self.jk_std = np.std(self.jk_data)

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

	def get_variance(self):
		"""
		Returns the variance of the jack knife method.
		"""
		self.var = 0
		for i in xrange(N):
			self.var += (jackknifed_dataset[i] - mean_biased)**2
		self.var = (N - 1) / float(N) * var
		return self.var

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


class Autocorrelation:
	"""
	Class for performing an autocorrelation analysis.
	"""
	def __init__(self, data):
		"""
		Args:
			data 					(numpy array): 	dataset to get autocorrelation for
		Returns:
			Object containing the autocorrelation values
		"""
		avg_data = np.average(data)
		self.N = len(data)
		self.data = data
		C0 = np.var(data)
		C = np.zeros(self.N/2)
		for h in xrange((self.N)/2):
			for i in xrange(0, self.N - h):
				C[h] += (data[i] - avg_data)*(data[i+h] - avg_data)
			C[h] /= (self.N - h)
		self.R = C / C0

	def __call__(self):
		"""
		Returns the auto-correlation
		"""
		return self.R

	def plot_autocorrelation(self, title, filename, lims = 1):
		"""
		Plots the autocorrelation.
		"""

		fig = plt.figure(dpi=300)
		ax = fig.add_subplot(111)
		ax.plot(range(self.N/2),self.R,color="r",label="Autocorrelation")

		# Giovanni function
		# ac2 = autocorrelation(self.data)
		# ax.plot(range(self.N/2),self.R,color="r",alpha=0.7,label="Autocorrelation")
		# ax.plot(range(self.N/2),ac2,color="b",alpha=0.7,label="Autocorrelation, numpy")

		ax.set_ylim(-lims,lims)
		ax.set_xlim(0,self.N/2)
		ax.set_xlabel(r"Lag $h$")
		ax.set_ylabel(r"$R = \frac{C_h}{C_0}$")
		ax.set_title(title)
		ax.grid(True)
		ax.legend()
		fig.savefig("autocorrelation_%s.png" % filename)

def autocorrelation(data):
    acf = np.zeros(len(data)/2)
    for k in range(0, len(data)/2):
        acf[k] = np.corrcoef(np.array([data[0:len(data)-k],data[k:len(data)]]))[0,1]
    return acf

def main():
	# Data to load and analyse
	data = np.loadtxt("data/prodRun1.dat",skiprows=8)
	
	# Histogram bins
	N_bins = 20
	
	# Bootstrapping
	N_bootstraps = int(1e4)
	bs = Bootstrap(data, N_bootstraps)
	bs_data = bs()

	# Jackknifing
	jk = Jackknife(data)
	jk_data = jk()

	# Autocorrelation
	ac = Autocorrelation(data)
	ac.plot_autocorrelation(r"Autocorrelation for $\beta = 6.0$", "beta6_0")

	print bs
	print jk

	plt.show()

if __name__ == '__main__':
	main()
