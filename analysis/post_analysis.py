import numpy as np, matplotlib.pyplot as plt, sys, os, re

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a # fermi

class DataReader:
	"""
	Small class for reading post analysis data
	"""
	def __init__(self,post_analysis_folder):
		self.post_analysis_folder = post_analysis_folder

		# Retrieves batch folders
		batch_folders = self._get_folder_content(post_analysis_folder)

		self.data_batches = {}
		for batch in batch_folders:
			batch_folder = os.path.join(post_analysis_folder,batch)
			batch_files = [os.path.join(batch_folder,i) for i in self._get_folder_content(batch_folder)]

			# Observable data list
			observable_data = []

			for bfile in batch_files:
				# Retrieves meta data from header
				meta_data = ""
				with open(bfile) as f:
					header_content = f.readline().split(" ")[1:]
					meta_data = [h.split("\n")[0] for h in header_content]

				# Retrieves observable from meta_data
				observable = meta_data[1] 

				observable_data.append({"observable": observable, "data": np.loadtxt(bfile)})

			# Loads data and stores it in a dictionary for all the data
			self.data_batches[batch] = observable_data

		self._reorganize_data()

	def _reorganize_data(self):
		# Reorganizes the data into beta-values and observables sorting
		self.data_observables = {}

		flow_time_retrieved = False

		# Sets up new dictionaries
		for batch_key in self.data_batches.keys():
			# print batch_key
			for obs_data in self.data_batches[batch_key]:
				self.data_observables[obs_data["observable"]] = {}

		# Places data into dictionaries
		for batch_key in self.data_batches.keys():
			for obs_data in self.data_batches[batch_key]:
				# Retrieves flow-time once
				if not flow_time_retrieved:
					self.flow_time = obs_data["data"][:,0]
					flow_time_retrieved = True

				self.data_observables[obs_data["observable"]][batch_key] = obs_data["data"][:,1:]

	def write_batch_to_single_file(self):
		file_folder = os.path.join(self.post_analysis_folder + "_universal_output")
		if not os.path.isdir(file_folder):
			os.mkdir(file_folder)
			print "> mkdir %s" % file_folder

		for batch_name in self.data_batches.keys():
			# Temporary data output for gathering data from arrays into a simple list
			_data_output = []
			header_output = ["t"]
			for obs_data in self.data_batches[batch_name]:
				# print obs_data["data"]
				_data_output.append(obs_data["data"][:,1:])
				# print obs_data
				header_output.append(obs_data["observable"])
				header_output.append(obs_data["observable"]+"_err")

			# Repopulating an array to write to file using numpy
			data_output = []
			for i in xrange(len(self.flow_time)):
				# Gathers the observables into a single list
				data_list = [self.flow_time[i]]
				for obs_data in _data_output: # Iterates over all the observables 
					for obs_data_val in obs_data[i]: # Iterates over value and error
						data_list.append(obs_data_val) # Appends all observables in the order of the header output

				data_output.append(data_list)			

			# Writing array to file
			file_path = os.path.join(file_folder,batch_name) + ".txt"
			np.savetxt(file_path,np.asarray(data_output),fmt="%.18f",header=" ".join(header_output))
			print "Batch '%s' data written to file." % batch_name

	@staticmethod
	def _get_folder_content(folder):
		if not os.path.isdir(folder):
			raise IOError("No folder by the name %s found." % folder)
		else:
			return os.listdir(folder)

class PostAnalysis:
	"""
	Post analysis base class
	"""
	observable_name = "Observable"
	observable_name_compact = "obs"
	x_label = r""
	y_label = r""
	dpi=None

	def __init__(self,data,flow_time,output_folder="../figures/post_analysis"):
		self.data = data
		self.flow_time = flow_time

		# Creates output folder for post analysis figures
		self.output_folder = output_folder
		if not os.path.isdir(self.output_folder):
			os.mkdir(self.output_folder)
			print "> mkdir %s" % self.output_folder

		# Creates colors to use
		self.colors = {}
		for color, batch_name in zip(["#5cbde0","#6fb718","#bc232e","#8519b7"],sorted(self.data.keys())): # blue, green, red purple
			self.colors[batch_name] = color

	@staticmethod
	def get_float(x):
		return float(".".join([i for i in re.findall('(\d+)',x)]))

	def _plot_core(self,plot_values,x_limits=False,y_limits=False):
		print "Plotting %s for betas %s together" % (self.analysis_name_compact," ".join(self.data.keys()))

		fig = plt.figure(dpi=self.dpi)
		ax = fig.add_subplot(111)

		for value in sorted(plot_values,key=lambda x:x["batch_name"]):
			x = value["x"]
			y = value["y"]
			y_err = value["y_err"]
			ax.plot(x,y,"-",label=value["label"],color=value["color"])
			ax.fill_between(x, y-y_err, y+y_err,alpha=0.5, edgecolor='', facecolor=value["color"])

		ax.grid(True)
		ax.set_title(r"%s %s" % (self.analysis_name, self.formula))
		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)

		if x_limits != False:
			ax.set_xlim(x_limits)
		if y_limits != False:
			ax.set_ylim(y_limits)

		ax.legend()

		fname = os.path.join(self.output_folder,"post_analysis_%s.png" % self.analysis_name_compact)
		plt.savefig(fname)
		print "Figure saved in %s" % fname

class TopSusPostAnalysis(PostAnalysis):
	analysis_name = "Topological Susceptibility"
	analysis_name_compact = "topsus"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"
	formula = r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$"

	def plot(self):
		plot_values = []
		for batch_name in sorted(self.data.keys()):
			values = {}
			values["batch_name"] = batch_name
			values["x"] = getLatticeSpacing(self.get_float(batch_name))*np.sqrt(8*self.flow_time)
			values["y"] = self.data[batch_name][:,0]
			values["y_err"] = self.data[batch_name][:,1]
			values["label"] = r"$\beta=%2.2f$" % self.get_float(batch_name)
			values["color"] = self.colors[batch_name]
			plot_values.append(values)

		self._plot_core(plot_values)

class EnergyPostAnalysis(PostAnalysis):
	analysis_name = "Energy"
	analysis_name_compact = "energy"
	y_label = r"$t^2\langle E\rangle$"
	x_label = r"$\frac{t}{r_0^2}$"
	formula = r"$t^2\langle E\rangle = t^2\frac{1}{64V}F_{\mu\nu}^a{F^{\mu\nu}}^a$"
	r0 = 0.5

	def plot(self):
		plot_values = []
		for batch_name in sorted(self.data.keys()):
			values = {}
			values["batch_name"] = batch_name
			values["x"] = self.flow_time/self.r0**2*getLatticeSpacing(self.get_float(batch_name))**2
			values["y"] = -self.data[batch_name][:,0]*self.flow_time**2/64.0
			values["y_err"] = self.data[batch_name][:,1]*self.flow_time**2/64.0
			values["label"] = r"$\beta=%2.2f$" % self.get_float(batch_name)
			values["color"] = self.colors[batch_name]
			plot_values.append(values)

		# x_limits = [0,np.max(plot_values[-1]["x"])]
		# y_limits = [0,np.max(plot_values[-1]["y"])]

		self._plot_core(plot_values)#,x_limits=x_limits,y_limits=y_limits)


def main(args):
	"""
	Args should be post-analysis folder
	"""
	# Loads data from post analysis folder
	data = DataReader(args[0])
	print "Retrieving data from folder: %s" % args[0]

	# Rewrites all of the data to a single file for sharing with giovanni
	data.write_batch_to_single_file()

	# Plots topsus
	topsus_analysis = TopSusPostAnalysis(data.data_observables["topsus"],data.flow_time)
	topsus_analysis.plot()

	# Retrofits the topsus for continiuum limit

	# Plots energy
	energy_analysis = EnergyPostAnalysis(data.data_observables["energy"],data.flow_time)
	energy_analysis.plot()

	# Retrofits the energy for continiuum limit
	plt.show()

if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		args = sys.argv[1:]
	else:
		args = ["../output/post_analysis_data"]
	main(args)