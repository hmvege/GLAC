import numpy as np, matplotlib.pyplot as plt, sys, os

class TopSusAnalysis:
	def __init__(self, folder_name):
		# Columns should consists of MC data, where each column is a datapoint
		datafiles = os.listdir(folder_name)
		datafiles = sorted(datafiles,key=lambda i : float(i.split('.out')[0].split('t')[-1])) # Must sort datafiles, if not will generate an error at t=2.0
		self.data = np.asarray([np.loadtxt(folder_name + '/' + filename)[:,1] for filename in datafiles])
		assert np.sum(np.asarray(np.loadtxt(folder_name + '/' + datafiles[0])[:,1]) - self.data[0]) < 1e-16, "Datapoints do not match up"
		self.N_data_points, self.N_sample_points = self.data.shape # Data points are points along the x-axis. Sample points are measurments at a single data point
		self.x = range(self.N_data_points)

	def set_x(self, x):
		"""
		For setting the x-axis points
		"""
		self.x = x

	def set_function(self, F):
		self.F = F

	def set_bs_statistics(self, Stat):
		"""
		Function for setting the statistics to be used.
		"""
		self.Stat = Stat

	def run_bootstraps(self, N_BS):
		"""
		Function that generates the new bootstrap dataset
		"""
		# Setting up arrays
		data_bs_raw = np.zeros((self.N_data_points,N_BS))
		# Index list for bootstrap - such that each data point uses the same index list
		index_lists = np.random.randint(self.N_sample_points, size=(N_BS, self.N_sample_points))
		# Performing the bootstrap samples
		for i in xrange(self.N_data_points):
			for j in xrange(N_BS):
				data_bs_raw[i,j] = self.Stat((self.data[i][index_lists[j]]))
		function_data = self.F(data_bs_raw)
		# Performing average and stds of the correlator
		self.bs 	= np.mean(function_data,axis=1)
		self.bs_std	= np.std(function_data,axis=1)
		# Updating class variables
		self.N_BS = N_BS
		self.N_sample_points = N_BS

	def plot(self,filename):
		plt.errorbar(self.x, y=self.bs, yerr=self.bs_std, fmt=".", ecolor="r", color="0",markevery=5,errorevery=5)
		plt.grid(True)
		plt.xlabel(r'$t_{flow}$',fontsize=18)
		plt.ylabel(r'$\chi_t^{1/4}[GeV]$',fontsize=18)
		plt.title(r'$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$, %s bootstraps' % self.N_BS, fontsize=18, y=1.025)
		plt.savefig('%s.png' % filename,dpi=300)
		print '%s.png written' % filename
		# plt.show()

	def plot_default(self,filename):
		plt.errorbar(self.x, y=np.mean(self.F(self.data**2),axis=1), yerr=np.std(self.F(self.data**2),axis=1), fmt=".", ecolor="r", color="0",markevery=5,errorevery=5)
		plt.grid(True)
		plt.xlabel(r'$t_{flow}$',fontsize=18)
		plt.ylabel(r'$\chi_t^{1/4}[GeV]$',fontsize=18)
		plt.title(r'$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$, no bootstraps', fontsize=18, y=1.025)
		plt.savefig('%s.png' % filename,dpi=300)
		print '%s.png written' % filename

	def write_bs_data(self,fname):		
		dat = np.zeros((len(self.x),2))
		for i in xrange(len(self.x)):
			dat[i,0] = self.x[i]
			dat[i,1] = self.bs[i]
		np.savetxt(fname,dat,fmt=["%10g","%15.15g"])
		print "%s written" % fname

	def write_average_data(self, fname):
		dat = np.zeros((len(self.x),2))
		for i in xrange(len(self.x)):
			dat[i,0] = self.x[i]
			dat[i,1] = self.F(self.Stat(self.data[i]))
		np.savetxt(fname,dat,fmt=["%10g","%15.15g"])
		print "%s written" % fname		

hbarc = 0.19732697 #eV micro m
a     = 0.0907 #fm
V     = 32**3 * 64
const = hbarc/a/V**(1./4)
const = 0.05717046003979148			# Jack's constant
# const = const / 0.05717046003979148

def chi(Q_squared):
	# Q should be averaged
	return const*Q_squared**(1./4)

def stat(x):
	return np.mean(x**2)

def main():
	folder_name = "data/PerTop/"
	N_BS 		= 1000
	filename 	= "figures/chi_%sbootstraps" % N_BS
	analysis = TopSusAnalysis(folder_name)
	analysis.set_x(np.arange(0,10,0.01) + 0.01)
	analysis.set_function(chi)
	analysis.plot_default("chi_default_plot")
	plt.show()
	exit(1)
	analysis.set_bs_statistics(stat)
	analysis.run_bootstraps(N_BS)
	analysis.write_bs_data("data/chi_analysis.txt")
	analysis.write_average_data("data/chi_non_bootstrap.txt")
	analysis.plot('figures/' + filename)

if __name__ == '__main__':
	main()