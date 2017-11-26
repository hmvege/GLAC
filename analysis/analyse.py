from LQCDAnalyser import Bootstrap, Jackknife, Autocorrelation
import os, numpy as np, matplotlib.pyplot as plt, sys, pandas as pd

class GetDirectoryTree:
	def __init__(self,batch_name):
		self.flow_tree = {}
		self.obs_tree = {}
		self.output_folder = "output"
		self.observables_list = ["plaq","topc","energy"]

		# Checks that the output folder actually exist
		if not os.path.isdir(os.path.join("..",self.output_folder)):
			raise EnvironmentError("No folder name output at location %s" % os.path.join("..",self.output_folder))

		# Retrieves folders and subfolders
		self.batch_folder = os.path.join("..",self.output_folder,batch_name)

		# Gets the regular configuration observables
		self.observables_folders = False
		if os.path.isdir(os.path.join(self.batch_folder,"observables")):
			self.observables_folder = os.path.join(self.batch_folder,"observables")
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
				# for i in self.flow_tree[flow_obs]:
				# 	print i

	def getFlow(self,obs):
		if obs in self.flow_tree.keys():
			return self.flow_tree[obs]

	def getObs(self,obs):
		if obs in self.obs_tree.keys():
			return self.obs_tree[obs]

class GetFolderContents:
	"""
	Analyses all the files in a given folder.
	"""
	def __init__(self,folder,flow=False):
		read_meta_data = True
		N_rows_to_skip = 0
		self.meta_data = {} # Assumes meta data is the same for all files in folder
		self.data_x = []
		self.data_y = []
		# Ensures we handle the data as a folder
		if type(folder) != list:
			folder = [folder]

		# Goes through files in folder and reads the contents into a file
		for file in folder:
			# Gets the metadata
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
				self.data_x.append(x)
			else:
				y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
			
			self.data_y.append(y)

		# Puts the data retrieved into arrays
		if flow: self.data_x = np.asarray(self.data_x) # Only provided in flow
		self.data_y = np.asarray(self.data_y)

def plot_flow_plaquette(x,y,y_error,meta_data,batch_name,N_bs=False):
	plt.figure()
	plt.errorbar(x,y,yerr=y_error,fmt="o",label="Plaquette")
	plt.xlabel(r"$\sqrt{8t_{flow}}$")
	plt.ylabel(r"$P_{\mu\nu}$")
	plt.title(r"$N_{flow}=%2d$, $\beta=%.2f$" % (meta_data["NFlows"],meta_data["beta"]))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/flow_plaquette_%s.png" % batch_name
	else:
		fname = "../figures/flow_plaquette_%s_Nbs%d.png" % (batch_name,N_bs)
	plt.savefig(fname)
	print "Figure created in %s" % fname

def plot_flow_topc(x,y,y_error,meta_data,batch_name,N_bs=False):
	plt.figure()
	plt.errorbar(x,y,yerr=y_error,fmt="o",label="Topological Charge")
	plt.xlabel(r"$\sqrt{8t_{flow}}$")
	plt.ylabel(r"$P_{\mu\nu}$")
	plt.title(r"$N_{flow}=%2d$, $\beta=%.2f$" % (meta_data["NFlows"],meta_data["beta"]))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/flow_topc_%s.png" % batch_name
	else:
		fname = "../figures/flow_topc_%s_Nbs%d.png" % (batch_name,N_bs)
	plt.savefig(fname)
	print "Figure created in %s" % fname

def plot_flow_energy(x,y,y_error,meta_data,batch_name,N_bs=False):
	plt.figure()
	plt.errorbar(x,y,yerr=y_error,fmt="o",label="Energy")
	plt.xlabel(r"$\sqrt{8t_{flow}}$")
	plt.ylabel(r"$\langle E \rangle$")
	plt.title(r"$N_{flow}=%2d$, $\beta=%.2f$" % (meta_data["NFlows"],meta_data["beta"]))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/flow_energy_%s.png" % batch_name
	else:
		fname = "../figures/flow_energy_%s_Nbs%d.png" % (batch_name,N_bs)
	plt.savefig(fname)
	print "Figure created in %s" % fname

def plot_flow_topsus(x,y,y_error,meta_data,batch_name,N_bs=False):
	plt.figure()
	plt.errorbar(x,y,yerr=y_error,fmt="o",label="Topological susceptibility")
	plt.xlabel(r"$\sqrt{8t_{flow}}$")
	plt.ylabel(r"$\chi_t^{1/4}[GeV]$")
	if not N_bs:
		plt.title(r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$, $\beta=%.2f$" % (meta_data["beta"]))
	else:
		plt.title(r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$, $\beta=%.2f$, $N_{bs}=%d$" % (meta_data["beta"],N_bs))
	plt.grid(True)
	if not N_bs:
		fname = "../figures/flow_topsus_%s.png" % batch_name
	else:
		fname = "../figures/flow_topsus_%s_Nbs%d.png" % (batch_name,N_bs)
	plt.savefig(fname)
	print "Figure created in %s" % fname

class AnalyseFlowObservable:
	def __init__(self,data):
		self.y = data.data_y
		self.N_configurations, self.NFlows = self.y.shape
		self.non_bs_data = np.zeros(self.NFlows)
		self.non_bs_data_std = np.zeros(self.NFlows)
		self.bs_data = np.zeros(self.NFlows)
		self.bs_data_std = np.zeros(self.NFlows)
		
	def boot(self,N_bs,B_statistic = np.mean, F = lambda x : x):
		for i in range(self.NFlows):
			bs = Bootstrap(self.y[:,i],N_bs, bootstrap_statistics = B_statistic, F = F, non_bs_stats = lambda x: x**2)
			self.non_bs_data[i] = bs.avg_original
			self.non_bs_data_std[i] = bs.std_original
			self.bs_data[i] = bs.bs_avg
			self.bs_data_std[i] = bs.bs_std

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a

hbarc = 0.19732697 #eV micro m
a     = 0.0931404061721 #fm
V     = 8**3 * 8
const = hbarc/a/V**(1./4)
# const = 0.05717046003979148			# Jack's constant
# const = const / 0.05717046003979148

def chi(Q_squared):
	# Q should be averaged
	return const*Q_squared**(1./4)

def stat(x,axis=None):
	return np.mean(x**2,axis=axis)


def main(args):
	if not args:
		args = ['verboseTest']

	DList = GetDirectoryTree(args[0])
	N_BS = 200

	# Flow observables
	flow_plaq_file = DList.getFlow("plaq")
	flow_topc_file = DList.getFlow("topc")
	flow_energy_file = DList.getFlow("energy")

	flow_plaq_data = GetFolderContents(flow_plaq_file,flow=True)
	flow_topc_data = GetFolderContents(flow_topc_file,flow=True)
	flow_energy_data = GetFolderContents(flow_energy_file,flow=True)

	flow_time = flow_topc_data.data_x[0]

	flow_plaq_analysis = AnalyseFlowObservable(flow_plaq_data)
	flow_topc_analysis = AnalyseFlowObservable(flow_topc_data)
	flow_energy_analysis = AnalyseFlowObservable(flow_energy_data)
	flow_topsus_analysis = AnalyseFlowObservable(flow_energy_data)
	
	flow_plaq_analysis.boot(N_BS)
	flow_topc_analysis.boot(N_BS)
	flow_energy_analysis.boot(N_BS)
	flow_topsus_analysis.boot(N_BS,stat)

	# No bootstrap
	plot_flow_plaquette(flow_time,flow_plaq_analysis.non_bs_data,flow_plaq_analysis.non_bs_data_std,flow_topc_data.meta_data,args[0])
	plot_flow_topc(flow_time,flow_topc_analysis.non_bs_data,flow_topc_analysis.non_bs_data_std,flow_topc_data.meta_data,args[0])
	plot_flow_energy(flow_time,flow_energy_analysis.non_bs_data,flow_energy_analysis.non_bs_data_std,flow_energy_data.meta_data,args[0])
	plot_flow_topsus(flow_time,flow_topsus_analysis.non_bs_data,flow_topsus_analysis.non_bs_data_std,flow_energy_data.meta_data,args[0])
	
	# With bootstrap
	plot_flow_plaquette(flow_time,flow_plaq_analysis.bs_data,flow_plaq_analysis.bs_data_std,flow_topc_data.meta_data,args[0],N_bs=N_BS)
	plot_flow_topc(flow_time,flow_topc_analysis.bs_data,flow_topc_analysis.bs_data_std,flow_topc_data.meta_data,args[0],N_bs=N_BS)
	plot_flow_energy(flow_time,flow_energy_analysis.bs_data,flow_energy_analysis.bs_data_std,flow_energy_data.meta_data,args[0],N_bs=N_BS)
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
