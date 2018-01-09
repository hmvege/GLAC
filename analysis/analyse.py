from LQCDAnalyser import Bootstrap, Jackknife, Autocorrelation
import os, numpy as np, matplotlib.pyplot as plt, sys, pandas as pd

class GetDirectoryTree:
	def __init__(self,batch_name,output_folder="output",dryrun=False):
		self.flow_tree = {}
		self.obs_tree = {}
		self.CURRENT_FOLDER = os.getcwd()
		self.output_folder = output_folder
		self.observables_list = ["plaq","topc","energy"]
		self.dryrun = dryrun

		# Checks that the output folder actually exist
		if not os.path.isdir(os.path.join("..",self.output_folder)):
			raise EnvironmentError("No folder name output at location %s" % os.path.join("..",self.output_folder))
		# Retrieves folders and subfolders
		self.batch_folder = os.path.join("..",self.output_folder,batch_name)

		# Gets the regular configuration observables
		self.observables_folders = False
		obs_path = os.path.join(self.batch_folder,"observables")
		if os.path.isdir(obs_path) and len(os.listdir(obs_path)) != 0:
			self.observables_folder = obs_path
			for obs,file_name in zip(self.observables_list,os.listdir(self.observables_folder)):
				obs_path = os.path.join(self.observables_folder,file_name)
				if os.path.isfile(obs_path):
					self.obs_tree[obs] = obs_path

		# Gets paths to flow observable
		# Checks that there exists a flow observable folder
		if os.path.isdir(os.path.join(self.batch_folder,"flow_observables")):
			# Creates the flow observables path
			flow_path = os.path.join(self.batch_folder,"flow_observables")
			# Goes through the flow observables
			for flow_obs in (self.observables_list):
				# Creates flow observables path
				obs_path = os.path.join(flow_path,flow_obs)
				# Checks if the flow observable path exists
				if os.path.isdir(obs_path):
					# Finds and sets the observable file paths
					flow_obs_dir_list = []
					for obs_file in os.listdir(obs_path):
						flow_obs_dir_list.append(os.path.join(obs_path,obs_file))
					self.flow_tree[flow_obs] = flow_obs_dir_list
		print "Directory tree built."

		# Creates figures folder
		self.figures_path = os.path.join("..","figures",batch_name)
		if not os.path.isdir(self.figures_path):
			if self.dryrun:
				print '> mkdir %s' % self.figures_path
			else:
				os.mkdir(self.figures_path)

	def getFlow(self,obs):
		"""
		Retrieves flow observable files.
		"""
		if obs in self.flow_tree.keys():
			return self.flow_tree[obs]
		else:
			raise Warning("Flow observable \"%s\" was not found in possible observables: %s" % (obs,", ".join(self.flow_tree.keys())))

	def getObs(self,obs):
		"""
		Retrieves observable files.
		"""
		if obs in self.obs_tree.keys():
			return self.obs_tree[obs]
		else:
			raise Warning("Observable \"%s\" was not found in possible observables: %s" % (obs,", ".join(self.flow_tree.keys())))

	def __str__(self):
		"""
		Prints the folder structre
		"""
		return_string = "Folder structure:"
		return_string += "\n{0:<s}".format(self.batch_folder)
		return_string += "\n{0:<s}/{1:<s}".format(self.batch_folder,"observables")
		if self.observables_folders:
			for obs,file_name in zip(self.observables_list,os.listdir(self.observables_folder)):
				return_string += "\n  {0:<s}".format(os.path.join(self.observables_folder,file_name))
		flow_path = os.path.join(self.batch_folder,"flow_observables")
		if os.path.isdir(flow_path):
			return_string += "\n  {0:<s}".format(flow_path)
			for flow_obs in (self.observables_list):
				obs_path = os.path.join(flow_path,flow_obs)
				return_string += "\n    {0:<s}".format(obs_path)
				for obs_file in os.listdir(obs_path):
					return_string += "\n      {0:<s}".format(os.path.join(obs_path,obs_file))
		return return_string

class GetFolderContents:
	"""
	Retrieves folder contents and acts as a container for data and meta-data.
	"""
	def __init__(self,folder,flow=False):
		if folder == None:
			print "No observables found in folder: %s" % folder
		else:
			read_meta_data = True
			retrieved_flow_time = False
			retrieved_indexing = False
			N_rows_to_skip = 0
			self.meta_data = {} # Assumes meta data is the same for all files in folder
			self.data_y = []
			self.data_x = False
			if flow:
				self.data_flow_time = False	
			# Ensures we handle the data as a folder
			if type(folder) != list:
				folder = [folder]

			# Goes through files in folder and reads the contents into a file
			for file in folder:
				# Gets the metadata
				print "    Reading file: %s" % file
				with open(file) as f:
					while read_meta_data:
						line = f.readline().split(" ")
						if line[0].isalpha():
							self.meta_data[str(line[0])] = float(line[-1])
							N_rows_to_skip += 1
						else:
							read_meta_data = False

				# Loads the data and places it in a list
				if flow:
					x, _x, y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
					if not retrieved_flow_time:
						self.data_flow_time = _x
				else:
					x, y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
				
				if not retrieved_indexing:
					self.data_x = x
				self.data_y.append(y)

			self.data_y = np.asarray(self.data_y)

			print "Data retrieved."


def plot_flow_plaquette(x,y,y_error,meta_data,batch_name,N_bs=False,dryrun=False):
	plt.figure()
	plt.errorbar(getLatticeSpacing(meta_data["beta"])*np.sqrt(8*x),y,yerr=y_error,fmt=".",color="0",ecolor="r",label="Plaquette",markevery=5,errorevery=5)
	plt.xlabel(r"$a\sqrt{8t_{flow}}$")
	plt.ylabel(r"$P_{\mu\nu}$")
	if not N_bs:
		plt.title(r"Plaquette $N_{flow}=%2d$, $\beta=%.2f$" % (meta_data["NFlows"],meta_data["beta"]))
	else:
		plt.title(r"Plaquette $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (meta_data["NFlows"],meta_data["beta"],N_bs))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/%s/flow_plaquette_%s.png" % (batch_name,batch_name)
	else:
		fname = "../figures/%s/flow_plaquette_%s_Nbs%d.png" % (batch_name,batch_name,N_bs)
	if not dryrun: plt.savefig(fname)
	print "Figure created in %s" % fname

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

class AnalyseFlowObservable:
	def __init__(self,data):
		self.y = data.data_y
		self.N_configurations, self.NFlows = self.y.shape
		# Non-bootstrapped data
		self.non_bs_data = np.zeros(self.NFlows)
		self.non_bs_data_std = np.zeros(self.NFlows)
		# Bootstrap data
		self.bs_data = np.zeros(self.NFlows)
		self.bs_data_std = np.zeros(self.NFlows)
		# Jackknifed data
		self.jk_data = np.zeros(self.NFlows)
		self.jk_data_std = np.zeros(self.NFlows)
		# Autocorrelation data
		self.autocorrelations = np.zeros((self.NFlows,self.N_configurations))

	def boot(self,N_bs,B_statistic = np.mean, F = lambda x : x, non_bs_stats = lambda x : x):
		for i in xrange(self.NFlows):
			bs = Bootstrap(self.y[:,i],N_bs, bootstrap_statistics = B_statistic, F = F, non_bs_stats = non_bs_stats)
			self.non_bs_data[i] = bs.avg_original
			self.non_bs_data_std[i] = bs.std_original
			self.bs_data[i] = bs.bs_avg
			self.bs_data_std[i] = bs.bs_std

	def jackknife(self,statistics = np.mean, F = lambda x : x):
		for i in xrange(self.NFlows):
			jk = Jackknife(self.y[i])
			self.jk_data[i] = jk.jk_avg
			self.jk_data_std[i] = jk.jk_std

	def autocorrelate(self):
		for i in xrange(self.NFlows):
			ac = Autocorrelation(self.y[i])
			self.autocorrelations[i] = ac()

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a

hbarc = 0.19732697 #eV micro m

# # Beta 6.0
# a     = getLatticeSpacing(6.0)
# V6_0     = 24**3 * 48
# const = hbarc/a/V6_0**(1./4)

# Beta 6.1
a     = getLatticeSpacing(6.1)
V6_1     = 28**3 * 56
const = hbarc/a/V6_1**(1./4)

def chi(Q_squared):
	# Q should be averaged
	return const*Q_squared**(1./4)

def stat(x,axis=None):
	return np.mean(x**2,axis=axis)

class Analyse(object):
	observable_name = "Missing_Observable_Name"
	x_label = "Missing x-label"
	y_label = "Missing y-label"
	mark_interval = 5
	error_mark_interval = 5

	def __init__(self,files,observable,batch_name,flow=True,data=None,dryrun=False):
		# Sets up global constants
		self.flow = flow
		self.batch_name = batch_name
		self.N_bs = None
		self.dryrun = dryrun

		print "Running analysis on run %s" % observable

		# Enables possibility of providing data
		if not data:
			self.data = GetFolderContents(files,flow=self.flow)
		else:
			self.data = data

		# Retrieves data
		if flow:
			self.x = self.data.data_x*self.data.meta_data["FlowEpsilon"]
		else:
			self.x = self.data.data_x

		# Creates an analyse object
		if flow:
			self.analysis = AnalyseFlowObservable(self.data)

		# Gets the lattice spacing
		self.beta = self.data.meta_data["beta"]
		self.a = getLatticeSpacing(self.beta)
		self.r = 0.5 # Sommer Parameters

	def boot(self,N_bs):
		self.N_bs = N_bs
		self.analysis.boot(N_bs)
		self.y = self.analysis.bs_data
		self.y_std = self.analysis.bs_data_std
		self.unanalyzed_y = self.analysis.non_bs_data
		self.unanalyzed_y_std = self.analysis.non_bs_data_std

	def jackknife(self):
		None

	def autocorrelation(self):
		None

	def plot(self,plot_bs=True):
		# Checks that the flow has been performed.
		if self.N_bs == None and plot_bs:
			raise ValueError("Flow has not been performed yet.")

		# Retrieves relevant data
		x = self.a * np.sqrt(8*self.x)
		if plot_bs:
			y = self.y
			y_std = self.y
		else:
			y = self.unanalyzed_y
			y_std = self.unanalyzed_y_std

		# Plotting commands
		plt.figure()
		plt.errorbar(x,y,yerr=y_std,fmt=".",color="0",ecolor="r",label=self.observable_name,markevery=self.mark_interval,errorevery=self.error_mark_interval)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label)
		plt.grid(True)
		if plot_bs:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta,self.N_bs)
			fname = "../figures/{0:<s}/flow_{2:<s}_{0:<s}_Nbs{1:<d}.png".format(self.batch_name,self.N_bs,"".join(self.observable_name.lower().split(" ")))
		else:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name, self.data.meta_data["NFlows"],self.beta)
			fname = "../figures/{0:<s}/flow_{1:<s}_{0:<s}.png".format(self.batch_name,"".join(self.observable_name.lower().split(" ")))
		plt.title(title_string)
		if not self.dryrun: plt.savefig(fname)
		print "Figure created in %s" % fname

class AnalysePlaquette(Analyse):
	observable_name = "Plaquette"
	x_label = r"$a\sqrt{8t_{flow}}$"
	y_label = r"$P_{\mu\nu}$"

	def __init__(self,files,observable,batch_name,flow=True,data=None,dryrun=False):
		super(AnalysePlaquette,self).__init__(files,observable,batch_name,flow=True,data=None,dryrun=False)

class AnalyseTopologicalCharge(Analyse):
	observable_name = "Topological Charge"
	x_label = r"$a\sqrt{8t_{flow}}$"
	y_label = r"$Q = \sum_x \frac{1}{32\pi^2}\epsilon_{\mu\nu\rho\sigma}Tr\{G^{clov}_{\mu\nu}G^{clov}_{\rho\sigma}\}$"

class AnalyseEnergy(Analyse):
	observable_name = "Energy"
	x_label = r"$t$"
	y_label = r"$\langle E \rangle t^2$"

	def plot(self,plot_bs=True):
		# Checks that the flow has been performed.
		if self.N_bs == None:
			raise ValueError("Flow has not been performed yet.")

		# Retrieves relevant data
		x = self.x/0.5**2
		if plot_bs:
			y = self.analysis.bs_data
			y_std = self.analysis.bs_data_std
		else:
			y = self.analysis.non_bs_data
			y_std = self.analysis.non_bs_data_std

		# Plotting commands
		plt.figure()
		plt.errorbar(x,y,yerr=y_std,fmt=".",color="0",ecolor="r",label=self.observable_name,markevery=self.mark_interval,errorevery=self.error_mark_interval)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label)
		plt.grid(True)
		if plot_bs:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$, $N_{bs}=%d$" % (self.observable_name,self.data.meta_data["NFlows"],self.beta,self.N_bs)
			fname = "../figures/{0:<s}/flow_{2:<s}_{0:<s}_Nbs{1:<d}.png".format(self.batch_name,self.N_bs,"".join(self.observable_name.lower().split(" ")))
		else:
			title_string = r"%s $N_{flow}=%2d$, $\beta=%.2f$" % (self.observable_name, self.data.meta_data["NFlows"],self.beta)
			fname = "../figures/{0:<s}/flow_{1:<s}_{0:<s}.png".format(self.batch_name,"".join(self.observable_name.lower().split(" ")))
		plt.title(title_string)
		if not self.dryrun: plt.savefig(fname)
		print "Figure created in %s" % fname

class AnalyseTopologicalSusceptibility(Analyse):
	observable_name = "Topological Susceptibility"
	x_label = r"$a\sqrt{8t_{flow}}$"
	y_label = r"$\chi_t^{1/4}[GeV]$"

	def plot(self,plot_bs=True):
		None

def main(args):
	if not args:
		args = ['prodRunBeta6_1','plaq','topc','energy','topsus']
		# args = ['prodRunBeta6_0','plaq','topc','energy','topsus']

	DList = GetDirectoryTree(args[0])
	N_bs = 200
	dryrun = False

	# Analyses plaquette data if present in arguments
	if 'plaq' in args:
		plaq_analysis = AnalysePlaquette(DList.getFlow("plaq"), "plaq", args[0], flow=True, dryrun = dryrun)
		plaq_analysis.boot(N_bs)
		plaq_analysis.plot()
		plaq_analysis.plot(plot_bs=False)

	sys.exit("EXITING: Missing complete plotter functions.")

	if 'topc' in args:
		topc_analysis = AnalyseTopologicalCharge(DList.getFlow("topc"), "topc", args[0], flow=True, dryrun = dryrun)
		topc_analysis.boot(N_bs)
		topc_analysis.plot()
		topc_analysis.plot(plot_bs=False)

		if 'topsus' in args:
			topsus_analysis = AnalyseTopologicalSusceptibility(DList.getFlow("topc"), "topsus", args[0], flow=True, dryrun = dryrun, data=topc_analysis.data)
			topsus_analysis.boot(N_bs)
			topsus_analysis.plot()
			topsus_analysis.plot(plot_bs=False)

	if 'energy' in args:
		energy_analysis = AnalyseEnergy(DList.getFlow("energy"), "energy", args[0], flow=True, dryrun = dryrun)
		energy_analysis.boot(N_bs)
		energy_analysis.plot()
		energy_analysis.plot(plot_bs=False)

	plt.show()
	# print flow_plaq_data.meta_data
	exit(1)

	# Flow observables
	flow_plaq_file = DList.getFlow("plaq")
	flow_topc_file = DList.getFlow("topc")
	flow_energy_file = DList.getFlow("energy") # Multiply by -1/4 on each clover, or, said in another way multiply by 1/16?

	flow_plaq_data = GetFolderContents(flow_plaq_file,flow=True)
	flow_topc_data = GetFolderContents(flow_topc_file,flow=True)
	flow_energy_data = GetFolderContents(flow_energy_file,flow=True)

	flow_time = flow_topc_data.data_x*0.01 # Scales to be in flow-time

	flow_plaq_analysis = AnalyseFlowObservable(flow_plaq_data)
	flow_topc_analysis = AnalyseFlowObservable(flow_topc_data)
	flow_energy_analysis = AnalyseFlowObservable(flow_energy_data)
	flow_topsus_analysis = AnalyseFlowObservable(flow_topc_data)
	
	flow_plaq_analysis.boot(N_BS)
	flow_topc_analysis.boot(N_BS)
	flow_energy_analysis.boot(N_BS)
	flow_topsus_analysis.boot(N_BS, B_statistic = stat, F = chi, non_bs_stats = lambda x : x**2)

	# Plaquette
	plot_flow_plaquette(flow_time,flow_plaq_analysis.non_bs_data,flow_plaq_analysis.non_bs_data_std,flow_topc_data.meta_data,args[0])
	plot_flow_plaquette(flow_time,flow_plaq_analysis.bs_data,flow_plaq_analysis.bs_data_std,flow_topc_data.meta_data,args[0],N_bs=N_BS)

	# Topological Charge
	plot_flow_topc(flow_time,flow_topc_analysis.non_bs_data,flow_topc_analysis.non_bs_data_std,flow_topc_data.meta_data,args[0])
	plot_flow_topc(flow_time,flow_topc_analysis.bs_data,flow_topc_analysis.bs_data_std,flow_topc_data.meta_data,args[0],N_bs=N_BS)

	# Action/energy density
	plot_flow_energy(flow_time,flow_energy_analysis.non_bs_data,flow_energy_analysis.non_bs_data_std,flow_energy_data.meta_data,args[0])
	plot_flow_energy(flow_time,flow_energy_analysis.bs_data,flow_energy_analysis.bs_data_std,flow_energy_data.meta_data,args[0],N_bs=N_BS)

	# Topological Susceptibility
	plot_flow_topsus(flow_time,flow_topsus_analysis.non_bs_data,flow_topsus_analysis.non_bs_data_std,flow_energy_data.meta_data,args[0])
	plot_flow_topsus(flow_time,flow_topsus_analysis.bs_data,flow_topsus_analysis.bs_data_std,flow_energy_data.meta_data,args[0],N_bs=N_BS)

	# Configuration observables
	plaq_file = DList.getObs("plaq")
	topc_file = DList.getObs("topc")
	energy_file = DList.getObs("energy")

	plaq_analysis = GetFolderContents(plaq_file)
	topc_analysis = GetFolderContents(topc_file)
	energy_analysis = GetFolderContents(energy_file)

if __name__ == '__main__':
	main(sys.argv[1:])
	plt.show()
