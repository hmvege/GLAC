from tools.folderreadingtools import GetDirectoryTree, GetFolderContents
from statistics.jackknife import Jackknife
from statistics.bootstrap import Bootstrap
from statistics.autocorrelation import Autocorrelation
import os, numpy as np, matplotlib.pyplot as plt, sys, pandas as pd

def plot_flow_topc(x,y,y_error,meta_data,batch_name,N_bs=False,dryrun=False):
	plt.figure()
	plt.errorbar(getLatticeSpacing(meta_data["beta"])*np.sqrt(8*x),y,yerr=y_error,fmt=".",color="0",ecolor="r",label="Topological Charge",markevery=5,errorevery=5)
	plt.xlabel(r"$a\sqrt{8t_{flow}}$")
	plt.ylabel(r"$Q = \sum_x \frac{1}{32\pi^2}\epsilon_{\mu\nu\rho\sigma}Tr\{G^{clov}_{\mu\nu}G^{clov}_{\rho\sigma}\}$")
	if not N_bs:
		plt.title(r"Topological Charge $N_{flow}=%2d$, $\beta=%.2f$" % (meta_data["NFlows"],meta_data["beta"]))
	else:
		plt.title(r"Topological Charge $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (meta_data["NFlows"],meta_data["beta"],N_bs))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/%s/flow_topc_%s.png" % (batch_name,batch_name)
	else:
		fname = "../figures/%s/flow_topc_%s_Nbs%d.png" % (batch_name,batch_name,N_bs)
	if not dryrun: plt.savefig(fname)
	print "Figure created in %s" % fname

def plot_flow_energy(x,y,y_error,meta_data,batch_name,N_bs=False,dryrun=False):
	plt.figure()
	plt.errorbar(x,-0.5*y*x**2,yerr=y_error,fmt=".",color="0",ecolor="r",label="Energy",markevery=5,errorevery=5)
	plt.xlabel(r"$t$")
	plt.ylabel(r"$\langle E \rangle t^2$")
	if not N_bs:
		plt.title(r"Energy Density $N_{flow}=%2d$, $\beta=%.2f$" % (meta_data["NFlows"],meta_data["beta"]))
	else:
		plt.title(r"Energy Density $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (meta_data["NFlows"],meta_data["beta"],N_bs))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/%s/flow_energy_%s.png" % (batch_name,batch_name)
	else:
		fname = "../figures/%s/flow_energy_%s_Nbs%d.png" % (batch_name,batch_name,N_bs)
	if not dryrun: plt.savefig(fname)
	print "Figure created in %s" % fname

def plot_flow_topsus(x,y,y_error,meta_data,batch_name,N_bs=False,dryrun=False):
	plt.figure()
	plt.errorbar(getLatticeSpacing(meta_data["beta"])*np.sqrt(8*x),y,yerr=y_error,fmt=".",color="0",ecolor="r",label="Topological susceptibility",markevery=5,errorevery=5)
	plt.xlabel(r"$a\sqrt{8t_{flow}}$")
	plt.ylabel(r"$\chi_t^{1/4}[GeV]$")
	if not N_bs:
		plt.title(r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$, $\beta=%.2f$" % (meta_data["beta"]))
	else:
		plt.title(r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$, $\beta=%.2f$, $N_{bs}=%d$" % (meta_data["beta"],N_bs))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/%s/flow_topsus_%s.png" % (batch_name,batch_name)
	else:
		fname = "../figures/%s/flow_topsus_%s_Nbs%d.png" % (batch_name,batch_name,N_bs)
	if not dryrun: plt.savefig(fname)
	print "Figure created in %s" % fname

class AnalyseFlow(object):
	observable_name = "Missing_Observable_Name"
	x_label = "Missing x-label"
	y_label = "Missing y-label"
	mark_interval = 5
	error_mark_interval = 5
	autocorrelations_limits = 1

	def __init__(self,files,observable,batch_name,data=None,dryrun=False,flow=True):
		# Sets up global constants
		self.batch_name = batch_name
		self.N_bs = None
		self.dryrun = dryrun
		self.flow = flow

		print "Running analysis on run %s" % observable

		# Enables possibility of providing data
		if data == None:
			self.data = GetFolderContents(files,flow=self.flow)
		else:
			self.data = data

		# Sets up variables
		self.y = self.data.data_y
		self.x = self.data.data_x
		self.N_configurations, self.NFlows = self.y.shape

		# Small error checking in retrieving number of flows
		if (int(self.data.meta_data["NFlows"]) + 1) != self.NFlows:
			raise ValueError("Number of flows %d does not match the number provided by metadata %d." % (self.NFlows,int(self.data.meta_data["NFlows"]) + 1))

		# Non-bootstrapped data
		self.unanalyzed_y = np.zeros(self.NFlows)
		self.unanalyzed_y_std = np.zeros(self.NFlows)

		# Bootstrap data
		self.bootstrap_performed = False
		self.bs_y = np.zeros(self.NFlows)
		self.bs_y_std = np.zeros(self.NFlows)

		# Jackknifed data
		self.jackknife_performed = False
		self.jk_y = np.zeros(self.NFlows)
		self.jk_y_std = np.zeros(self.NFlows)

		# Autocorrelation data
		self.autocorrelation_performed = False
		self.autocorrelations = np.zeros((self.NFlows,self.N_configurations / 2))

		# Gets the lattice spacing
		self.beta = self.data.meta_data["beta"]
		self.a = self.getLatticeSpacing(self.beta)
		self.r = 0.5 # Sommer Parameters

	@staticmethod
	def getLatticeSpacing(beta):
		# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
		r = 0.5
		bval = (beta - 6)
		a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
		return a

	def boot(self,N_bs,B_statistic = np.mean, F = lambda x : x, non_bs_stats = lambda x : x):
		self.N_bs = N_bs
		index_lists = np.random.randint(self.N_configurations, size=(self.N_bs, self.N_configurations))
		for i in xrange(self.NFlows):
			bs = Bootstrap(self.y[:,i],N_bs, bootstrap_statistics = B_statistic, F = F, non_bs_stats = non_bs_stats, index_lists = index_lists)
			self.bs_y[i] = bs.bs_avg
			self.bs_y_std[i] = bs.bs_std
			self.unanalyzed_y[i] = bs.avg_original
			self.unanalyzed_y_std[i] = bs.std_original

		# Sets performed flag to true
		self.bootstrap_performed = True

	def jackknife(self,data_statistics = np.mean, F = lambda x : x):
		for i in xrange(self.NFlows):
			jk = Jackknife(self.y[:,i],F = F, data_statistics = data_statistics)
			self.jk_y[i] = jk.jk_avg
			self.jk_y_std[i] = jk.jk_std

		# Sets performed flag to true
		self.jackknife_performed = True

	def autocorrelation(self):
		for i in xrange(self.NFlows):
			ac = Autocorrelation(self.y[:,i])
			self.autocorrelations[i] = ac()
			# # Small progressbar
			# sys.stdout.write("\r%3.1f%% done" % (100*float(i)/float(self.NFlows)))
			# sys.stdout.flush()

		# Sets performed flag to true
		self.autocorrelation_performed = True

	def plot_jackknife(self, x = None, y_function = lambda x : x):
		# Checks that jacknifing has been performed.
		if not self.jackknife_performed:
			raise ValueError("Jackknifing has not been performed yet.")

		# Sets up the x axis array to be plotted
		if not x:
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x*self.data.meta_data["FlowEpsilon"])

		# Sets up the title and filename strings
		title_string = r"Jacknife of %s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta)
		fname = "../figures/{0:<s}/flow_{1:<s}_{0:<s}_jackknife.png".format(self.batch_name,"".join(self.observable_name.lower().split(" ")))

		# Plots the jackknifed data
		plt.figure()
		plt.errorbar(x,self.jk_y,yerr=self.jk_y_std,fmt=".",color="0",ecolor="r",label=self.observable_name,markevery=self.mark_interval,errorevery=self.error_mark_interval)
		plt.xlabel(self.x_label)
		plt.xlabel(self.y_label)
		plt.title(title_string)
		if not self.dryrun: 
			plt.savefig(fname,dpi=300)
		print "Figure created in %s" % fname

	def plot_autocorrelation(self):
		# Checks that autocorrelations has been performed.
		if not self.autocorrelation_performed:
			raise ValueError("Autocorrelation has not been performed yet.")

		# sets up the autocorrelation
		N_autocorr = self.N_configurations / 2

		# Finds the maximum value at each MC time and sets up the y array
		x = range(N_autocorr)
		y = np.zeros(N_autocorr)
		for i in xrange(N_autocorr):
			y[i] = np.max(self.autocorrelations[:,i])

		# Sets up the title and filename strings
		title_string = r"Autocorrelation of %s $N_{flow}=%2d$, $\beta=%.2f$, $N_{cfg}=%2d$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta,self.N_configurations)
		fname = "../figures/{0:<s}/flow_{1:<s}_{0:<s}_autocorrelation.png".format(self.batch_name,"".join(self.observable_name.lower().split(" ")))

		# Plots the autocorrelations
		fig = plt.figure(dpi=300)
		ax = fig.add_subplot(111)
		ax.plot(x,np.abs(y),color="0",label=self.observable_name)
		ax.set_ylim(-self.autocorrelations_limits,self.autocorrelations_limits)
		ax.set_xlim(0,N_autocorr)
		ax.set_xlabel(r"Lag $h$")
		ax.set_ylabel(r"$R = \frac{C_h}{C_0}$")
		ax.set_title(title_string)
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.2))
		ax.grid(True)
		ax.legend()
		if not self.dryrun: 
			fig.savefig(fname,dpi=300)
		print "Figure created in %s" % fname

	def plot_boot(self,plot_bs=True, x = None, y_function = lambda x : x):
		# Checks that the bootstrap has been performed.
		if not self.bootstrap_performed and plot_bs:
			raise ValueError("Bootstrap has not been performed yet.")

		# Retrieves relevant data and sets up the arrays to be plotted
		if not x:
			# Default x axis points is the flow time
			x = self.a * np.sqrt(8*self.x*self.data.meta_data["FlowEpsilon"])
		if plot_bs:
			y = self.bs_y
			y_std = self.bs_y_std
		else:
			y = self.unanalyzed_y
			y_std = self.unanalyzed_y_std

		# Sets up the title and filename strings
		if plot_bs:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta,self.N_bs)
			fname = "../figures/{0:<s}/flow_{2:<s}_{0:<s}_Nbs{1:<d}.png".format(self.batch_name,self.N_bs,"".join(self.observable_name.lower().split(" ")))
		else:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name, self.data.meta_data["NFlows"],self.beta)
			fname = "../figures/{0:<s}/flow_{1:<s}_{0:<s}.png".format(self.batch_name,"".join(self.observable_name.lower().split(" ")))

		# Plots either bootstrapped or regular stuff
		plt.figure()
		plt.errorbar(x,y,yerr=y_std,fmt=".",color="0",ecolor="r",label=self.observable_name,markevery=self.mark_interval,errorevery=self.error_mark_interval)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label)
		plt.grid(True)
		plt.title(title_string)
		if not self.dryrun:
			plt.savefig(fname,dpi=300)
		print "Figure created in %s" % fname

	def plot_original(self):
		"""
		Plots the default analysis, mean and std of the observable.
		"""
		self.plot_boot(plot_bs=False)

class AnalysePlaquette(AnalyseFlow):
	"""
	Plaquette analysis class.
	"""
	observable_name = "Plaquette"
	x_label = r"$a\sqrt{8t_{flow}}[fm]$"
	y_label = r"$P_{\mu\nu}$"

class AnalyseTopologicalCharge(AnalyseFlow):
	"""
	Topological charge analysis class. NOT TESTED
	"""
	observable_name = "Topological Charge"
	x_label = r"$a\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q = \sum_x \frac{1}{32\pi^2}\epsilon_{\mu\nu\rho\sigma}Tr\{G^{clov}_{\mu\nu}G^{clov}_{\rho\sigma}\}$"

class AnalyseEnergy(AnalyseFlow):
	"""
	Energy/action density analysis class. NOT TESTED
	"""
	observable_name = "Energy"
	x_label = r"$t/r_0^2$"
	y_label = r"$\langle E \rangle t^2$" # Energy is dimension 4, while t^2 is dimension invsere 4, or length/time which is inverse energy, see Peskin and Schroeder

class AnalyseTopologicalSusceptibility(AnalyseFlow):
	"""
	Topological susceptibility analysis class. NOT TESTED / NOT COMPLETED
	"""
	observable_name = "Topological Susceptibility"
	x_label = r"$a\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]"

	def __init__(self,files,observable,batch_name,data=None,dryrun=False):
		super(AnalyseTopologicalSusceptibility,self).__init__(files,observable,batch_name,data=None,dryrun=False)
		self.set_size()

	def set_size(self):
		self.a = getLatticeSpacing(meta_data["beta"])
		
		# Ugly, hardcoded lattice size setting
		if meta_data["beta"] == 6.0:
			lattice_size = 24**3*48
		else if meta_data["beta"] == 6.1:
			lattice_size = 28**3*56
		else if meta_data["beta"] == 6.2:
			lattice_size = 32**3*64
		else if meta_data["beta"] == 6.45:
			lattice_size = 48**3*96
		else:
			raise ValueError("Unrecognized beta value: %g" % meta_data["beta"])

		self.V = lattice_size
		self.hbarc = 0.19732697 #eV micro m
		self.const = self.hbarc/self.a/self.V**(1./4)

	@staticmethod
	def chi(Q_squared):
		# Q should be averaged
		return self.const*Q_squared**(1./4)

	@staticmethod
	def stat(x,axis=None):
		return np.mean(x**2,axis=axis)

def main(args):
	if not args:
		# args = ['prodRunBeta6_1','plaq','topc','energy','topsus']
		# args = ['prodRunBeta6_0','plaq','topc','energy','topsus']
		args = ['beta6_1','plaq','topc','energy','topsus']

	DirectoryList = GetDirectoryTree(args[0],output_folder="data")
	N_bs = 200
	dryrun = False
	# print DirectoryList

	# Analyses plaquette data if present in arguments
	# if 'plaq' in args:
	# 	plaq_analysis = AnalysePlaquette(DirectoryList.getFlow("plaq"), "plaq", args[0], dryrun = dryrun)
	# 	plaq_analysis.boot(N_bs)
	# 	plaq_analysis.jackknife()
	# 	plaq_analysis.autocorrelation()
	# 	plaq_analysis.plot_boot()
	# 	plaq_analysis.plot_original()
	# 	plaq_analysis.plot_jackknife()
	# 	plaq_analysis.plot_autocorrelation()

	if 'topc' in args:
		topc_analysis = AnalyseTopologicalCharge(DirectoryList.getFlow("topc"), "topc", args[0], dryrun = dryrun)
		# topc_analysis.boot(N_bs)
		# topc_analysis.jackknife()
		# topc_analysis.autocorrelation()
		# topc_analysis.plot_boot()
		# topc_analysis.plot_original()
		# topc_analysis.plot_jackknife()
		# topc_analysis.plot_autocorrelation()


		if 'topsus' in args:
			topsus_analysis = AnalyseTopologicalSusceptibility(DirectoryList.getFlow("topc"), "topsus", args[0], dryrun = dryrun, data=topc_analysis.data)
			topsus_analysis.boot(N_bs,B_statistic = topsus_analysis.stat, F = topsus_analysis.chi, non_bs_stats = topsus_analysis.stat)
			topsus_analysis.jackknife(data_statistics = lambda x : np.mean(x**2), F = topsus_analysis.chi)
			# topsus_analysis.autocorrelation() # Dosen't make sense to do the autocorrelation of the topoligical susceptibility since it is based on a data mean
			topsus_analysis.plot_boot()
			topsus_analysis.plot_original()
			topsus_analysis.plot_jackknife()
			topsus_analysis.plot_autocorrelation()

		sys.exit("EXITING: Missing complete plotter functions.")

	if 'energy' in args:
		energy_analysis = AnalyseEnergy(DirectoryList.getFlow("energy"), "energy", args[0], dryrun = dryrun)
		energy_analysis.boot(N_bs)
		energy_analysis.jackknife()
		energy_analysis.autocorrelation()
		energy_analysis.plot_boot()
		energy_analysis.plot_original()
		energy_analysis.plot_jackknife()
		energy_analysis.plot_autocorrelation()

	# plt.show()
	# print flow_plaq_data.meta_data
	exit(1)

if __name__ == '__main__':
	main(sys.argv[1:])
	plt.show()
