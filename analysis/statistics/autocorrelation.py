import numpy as np, matplotlib.pyplot as plt, sys, os

__all__ = ["Autocorrelation"]

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

	def plot_autocorrelation(self, title, filename, lims = 1,dryrun=False):
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
		if dryrun:
			fig.savefig("autocorrelation_%s.png" % filename)

def alternate_autocorrelation(data):
    acf = np.zeros(len(data)/2)
    for k in range(0, len(data)/2):
        acf[k] = np.corrcoef(np.array([data[0:len(data)-k],data[k:len(data)]]))[0,1]
    return acf

def main():
	# Data to load and analyse
	data = np.loadtxt("tests/plaq.dat",skiprows=8)
	
	# Histogram bins
	N_bins = 20
	
	# Autocorrelation
	ac = Autocorrelation(data)
	ac.plot_autocorrelation(r"Autocorrelation for $\beta = 6.1$", "beta6_1",dryrun=True)

	plt.show()

if __name__ == '__main__':
	main()