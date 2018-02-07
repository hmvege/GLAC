import numpy as np, matplotlib.pyplot as plt, sys, os, re, scipy.optimize as sciopt, copy
from tools.folderreadingtools import check_folder

# import tqdm
# from scipy.optimize import curve_fit

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a # fermi

class PostAnalysisDataReader:
	"""
	Small class for reading post analysis data
	"""
	def __init__(self,post_analysis_folder,verbose=False):
		self.post_analysis_folder = post_analysis_folder
		self.verbose = verbose

		# Dictionary variable to hold all the data sorted by batches
		self.data_batches = {}

		# Dictionaries to hold the raw bootstrapped/jackknifed data and autocorrelation data
		self.bs_data_raw = {}
		self.jk_data_raw = {}
		self.ac_corrections = {}

		# Different types of analysis
		self.analysis_types = ["unanalyzed","jackknife","bootstrap"]

		# Binary folder types available
		self.binary_folder_types = ["jackknife","bootstrap","autocorrelation"]

		# Variable to store if we have retrieved flow time or not
		self.retrieved_flow_time = False

		# Data batch name
		self.data_batch_name = os.path.split(self.post_analysis_folder)[-1]

		# Iterates over the different beta value folders
		for beta_folder in self._get_folder_content(self.post_analysis_folder):
			# Construct beta folder path
			beta_folder_path = os.path.join(self.post_analysis_folder,beta_folder)

			# # Stores the batch folder path
			# batch_folder = os.path.join(self.post_analysis_folder,beta_batch)

			# # Stores retrieves files in batch folder
			# batch_files = [os.path.join(batch_folder,i) for i in self._get_folder_content(batch_folder)]

			# Retrieves beta from folder name
			beta = float(beta_folder.strip("beta").replace("_","."))

			# Dictionary to store observable data in
			observable_data = {}

			# Loops over files and folders inside the beta folder
			for folder in self._get_folder_content(beta_folder_path):
				# Construct folder path
				folder_path = os.path.join(beta_folder_path,folder)

				# Retrieves data deepending on type
				if folder == "bootstrap":
					# Gets binary bs data
					self.bs_data_raw[beta] = self._get_bin_dict(folder_path)
				elif folder == "jackknife":
					# Gets binary jk data
					self.jk_data_raw[beta] = self._get_bin_dict(folder_path)
				elif folder == "autocorrelation":
					# Gets binary ac data
					self.ac_corrections[beta] = self._get_bin_dict(folder_path)
				else:
					# Gets observable name
					observable_name = os.path.splitext(folder)[0]
					
					# Retrieves the observable data
					observable_data[observable_name] = self._get_beta_observable_dict(folder_path,beta)

			# Stores batch data
			self.data_batches[beta] = copy.deepcopy(observable_data)

			# Frees memory
			del observable_data
		
		# Reorganizes data to more ease-of-use type of data set
		self._reorganize_data()

	def _get_beta_observable_dict(self,observable_file,beta):
		# Retrieves metadata
		# Make it so one can retrieve the key as meta_data[i] and then value as meta_data[i+1]
		meta_data = self._get_meta_data(observable_file)

		# Temporary methods for getting observable name and beta value, as this will be put into the meta data
		obs = os.path.split(os.path.splitext(observable_file)[0])[-1]

		# Dictionary to store all observable data in
		obs_data = {}

		# Loads data into temporary holder
		retrieved_data = np.loadtxt(observable_file)

		# Puts data into temporary holding facilities
		t 			= retrieved_data[:,0]
		y 			= retrieved_data[:,1]
		y_error 	= retrieved_data[:,2]
		bs_y 		= retrieved_data[:,3]
		bs_y_error 	= retrieved_data[:,4]
		jk_y 		= retrieved_data[:,5]
		jk_y_error 	= retrieved_data[:,6]


		# Stores data into dictionaries
		unanalyzed_data = {"y": y, "y_error": y_error}
		bs_data = {"y": bs_y, "y_error": bs_y_error}
		jk_data = {"y": jk_y, "y_error": jk_y_error}

		# Stores observable data
		obs_data["beta"] 		= copy.deepcopy(beta)
		obs_data["unanalyzed"] 	= copy.deepcopy(unanalyzed_data)
		obs_data["bootstrap"] 	= copy.deepcopy(bs_data)
		obs_data["jackknife"] 	= copy.deepcopy(jk_data)

		# Stores flow time in a seperate variable
		if not self.retrieved_flow_time:
			self.flow_time = copy.deepcopy(t)
			self.retrieved_flow_time = True

		if self.verbose:
			print "Data retrieved from %s" % observable_file

		# Frees memory
		del retrieved_data

		return obs_data

	def _get_bin_dict(self,folder):
		"""
		Gets binary data files
		"""
		observable_dict = {}

		# Retrieves the different bootstrapped observables
		for observable_file in self._get_folder_content(folder):
			# Gets observable name
			observable_name = os.path.splitext(observable_file)[0]

			# Populates dictionary
			observable_dict[observable_name] = np.load(os.path.join(folder,observable_file))

		return observable_dict

	@staticmethod
	def _get_meta_data(file):
		# Retrieves meta data from header or file
		meta_data = ""
		with open(file) as f:
			header_content = f.readline().split(" ")[1:]
			meta_data = [h.split("\n")[0] for h in header_content]
		return meta_data

	def _reorganize_data(self):
		# Reorganizes the data into beta-values and observables sorting
		self.data_observables = {}

		# Sets up new dictionaries by looping over batch names
		for beta in self.data_batches:
			# Loops over observable names
			for observable_name in self.data_batches[beta]:
				# Creates new sub-dictionary ordered by the observable name
				self.data_observables[observable_name] = {}

		# Places data into dictionaries
		for beta in self.data_batches:
			# Loops over the batch observable
			for observable_name in self.data_batches[beta]:
				# Stores the batch data in a sub-dictionary
				self.data_observables[observable_name][beta] = self.data_batches[beta][observable_name]


	def write_batch_to_single_file(self):
		"""
		Writes all unanalyzed data to a single file
		"""
		# Creates universal output folder
		file_folder = os.path.join(self.post_analysis_folder,"..","..","universal_output") # Should be in output
		print "="*100,"\nWriting data to a universal output format in %s" % file_folder
		check_folder(file_folder,False,True)

		# Variable for checking if we have retrieved the flow time or not
		retrieved_flow_time = False

		# Loops over the different possible analyses
		for analysis_type in self.analysis_types:
			# Loops over the batches
			for beta in self.data_batches:
				# Temporary data output for gathering data from arrays into a simple list
				data_output = []

				# Header for the resulting array
				header_output = ["t","sqrt8t"]

				# Retrieves t
				data_output.append(self.flow_time)

				# Retrieves sqrt(8t)
				data_output.append(np.sqrt(8*self.flow_time))

				for observable_name in self.data_batches[beta]:
					# Retrieves the unanalyzed observable data
					data_output.append(self.data_batches[beta][observable_name][analysis_type]["y"])
					data_output.append(self.data_batches[beta][observable_name][analysis_type]["y_error"])

					# Adds observable name and error to header
					header_output.append(observable_name)
					header_output.append(observable_name+"_err")

				# Writing array to file
				file_path = os.path.join(file_folder,str(beta).replace(".","_") + "_" + analysis_type) + ".txt"
				np.savetxt(file_path,np.asarray(data_output),fmt="%.18f",header=" ".join(header_output))
				print "Batch '%s' data written to file %s." % (beta,os.path.splitext(os.path.split(file_path)[-1])[0])

		# A single delimiter marking the end of universal output writing
		print "="*100

	@staticmethod
	def _get_folder_content(folder):
		if not os.path.isdir(folder):
			raise IOError("No folder by the name %s found." % folder)
		else:
			return [ f for f in os.listdir(folder) if not f.startswith(".") ]

def fit_line(x, bs_data_dict, fit_target, fit_interval, data_function_correction = lambda x : x, fit_function_modifier = lambda x : x, fit_method = "curve_fit", plot_fit_window = False):
	"""
	Function for creating a line fit using bootstrapped data.
	Args:
		x 					(numpy float array): x-axis values to fit along
		bs_data_dict  		 	   (dictionary): dictionary containing bootstrapped array "data", "observable" and "beta"
		fit_target						(float): y-axis value to fit against
		fit_interval					(float): y-axis interval in which we will linefit for
		[data_function_correction]	 (function): function for correcting data before analysis is performed
		[fit_function_modifier] 	 (function): function for modifying y-values
		[fit_method]					  (str): type of data fit to be performed. Default is scipy.optimize.curve_fit
		[plot_fit_window]				 (bool): plots the fit window to easier visualize
	Returns:
		fit_value 						(float): y value fitted at fit_target
		fit_value_error 				(float): error of y_0
	"""
	# Retrieves data
	print "@175 ",bs_data_dict.keys()
	observable = bs_data_dict["observable"]
	beta = bs_data_dict["beta"]
	data = bs_data_dict["data"]
	NFlows, NBoot = data.shape


	exit(1)
	# Applies the funciton correction to all bootstrap samples
	for iBoot in xrange(NBoot):
		data[:,iBoot] = data_function_correction(data[:,iBoot])

	# Sets the fitting method for the each bootstrap value
	fit_methods_list = ["curve_fit","polynomial"]
	if fit_method not in fit_methods_list:
		raise KeyError("%s not a possible fit method. Use %s" % (fit_method,", ".join(fit_methods_list)))

	# Takes data mean for finding fit interval
	data_mean = np.mean(bs_data_dict["data"],axis=1)
	data_std = np.std(bs_data_dict["data"],axis=1)

	# Finds fit interval indices
	start_index = np.argmin( np.abs(data_mean - (fit_target - fit_interval)) )
	end_index = np.argmin( np.abs(data_mean - (fit_target + fit_interval)) )
	N_fit_line_length = end_index - start_index

	# Plots the fit window if prompted to do so
	if plot_fit_window:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(x,data_mean,label=observable)
		ax.fill_between(x, data_mean - data_std, data_mean + data_std,alpha=0.5)
		ax.scatter(x[start_index],data_mean[start_index],label="Fit line start")
		ax.scatter(x[end_index],data_mean[end_index],label="Fit line end")
		ax.axhline(fit_target)
		ax.grid(True)
		ax.legend()
		ax.set_title("Fit window for %s" % observable)
		plt.show()
		plt.close(fig)

	# Array to store fitted values in
	fitted_values = np.zeros(NBoot)

	for iBoot in xrange(NBoot):
		# Fitting data
		if fit_method == "curve_fit":
			pol,polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x[start_index:end_index],data[start_index:end_index,iBoot])
		elif fit_method == "polynomial":
			pol,polcov = np.polyfit(x[start_index:end_index],data[start_index:end_index,iBoot],1,rcond=None,full=False,cov=True)
		else:
			raise KeyError("No fit method called %s found." % fit_method)

		# Retrieving polynomial values for retrofitting
		a = pol[0]
		b = pol[1]

		# Adds bootstrap fit to fitted values
		fitted_values[iBoot] = fit_function_modifier((fit_target - b) / a)

	# Takes mean and standard deviation of fitted values
	fit_value = np.mean(fitted_values)
	fit_value_error = np.std(fitted_values)

	return fit_value, fit_value_error


class PostAnalysis:
	"""
	Post analysis base class
	"""
	observable_name = "Observable"
	observable_name_compact = "obs"
	x_label = r""
	y_label = r""
	dpi=None
	size_labels = {	6.0  : r"$24^3 \times 48$",
					6.1  : r"$28^3 \times 56$",
					6.2  : r"$32^3 \times 64$",
					6.45 : r"$48^3 \times 96$"}

	# def __init__(self,data,flow_time,data_batch_name,base_output_folder="../figures/post_analysis"):
	def __init__(self,data,observable,base_output_folder="../figures/post_analysis"):
		# Retrieves relevant data values
		self.flow_time 			= data.flow_time
		self.unanalyzed_data 	= {beta:data.data_observables[observable][beta]["unanalyzed"] for beta in sorted(data.data_observables[observable].keys())} # Data should now be sorted by beta values
		self.bootstrap_data 	= {beta:data.data_observables[observable][beta]["bootstrap"] for beta in sorted(data.data_observables[observable].keys())}
		self.jackknife_data 	= {beta:data.data_observables[observable][beta]["jackknife"] for beta in sorted(data.data_observables[observable].keys())}
		self.bs_raw 			= data.bs_data_raw
		self.jk_raw 			= data.jk_data_raw
		self.ac_corrections		= data.ac_corrections

		# Creates base output folder for post analysis figures
		self.base_output_folder_path = base_output_folder
		check_folder(self.base_output_folder_path,dryrun=False,verbose=True)

		# Creates output folder
		self.output_folder_path = os.path.join(self.base_output_folder_path,data.data_batch_name)
		check_folder(self.output_folder_path,dryrun=False,verbose=True)

		# Creates colors to use
		self.colors = {}
		for color, beta in zip(["#5cbde0","#6fb718","#bc232e","#8519b7"],sorted(data.data_observables[observable].keys())): # blue, green, red purple
			self.colors[beta] = color

	@staticmethod
	def get_float(x):
		return float(".".join([i for i in re.findall('(\d+)',x)]))

	def _plot_core(self, x_limits = False, y_limits = False, fname_appendix = ""):
		print "Plotting %s for betas %s together" % (self.observable_name_compact,", ".join([str(b) for b in sorted(self.unanalyzed_data.keys())]))
		fig = plt.figure(dpi=self.dpi)
		ax = fig.add_subplot(111)

		for value in self.plot_values:
			x = value["x"]
			y = value["y"]
			y_err = value["y_err"]
			ax.plot(x,y,"-",label=value["label"],color=value["color"])
			ax.fill_between(x, y-y_err, y+y_err,alpha=0.5, edgecolor='', facecolor=value["color"])

		ax.grid(True)
		ax.set_title(r"%s %s" % (self.observable_name, self.formula))
		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)

		if x_limits != False:
			ax.set_xlim(x_limits)
		if y_limits != False:
			ax.set_ylim(y_limits)

		ax.legend(loc="lower right")

		fname = os.path.join(self.output_folder_path,"post_analysis_%s%s.png" % (self.observable_name_compact,fname_appendix))
		plt.savefig(fname)
		print "Figure saved in %s" % fname
		# plt.show()
		plt.close(fig)

class TopSusPostAnalysis(PostAnalysis):
	observable_name = "Topological Susceptibility"
	observable_name_compact = "topsus"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"
	formula = r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$"

	def __init__(self,*args,**kwargs):
		super(TopSusPostAnalysis,self).__init__(*args,**kwargs)

	def plot(self,analysis_type="bootstrap"):
		self.plot_values = []

		# Retrieving data deepending on analysis type we are choosing
		if analysis_type == "bootstrap":
			data = self.bootstrap_data
		elif analysis_type == "jackknife":
			data = self.jackknife_data
		elif analysis_type == "unanalyzed":
			data = self.unanalyzed_data
		else:
			raise KeyError("Analysis %s not recognized" % analysis_type)

		# Sorts data into a format specific for the plotting method
		for beta in sorted(data.keys()):
			values = {}
			values["beta"] = beta
			values["x"] = np.sqrt(8*self.flow_time)
			values["y"] = data[beta]["y"]
			values["y_err"] = data[beta]["y_error"]
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta],beta)
			values["color"] = self.colors[beta]
			self.plot_values.append(values)

		self._plot_core(fname_appendix = "_" + analysis_type)

class EnergyPostAnalysis(PostAnalysis):
	observable_name = "Energy"
	observable_name_compact = "energy"
	y_label = r"$t^2\langle E\rangle$"
	x_label = r"$\frac{t}{r_0^2}$"
	formula = r"$\langle E\rangle = -\frac{1}{64V}F_{\mu\nu}^a{F^a}^{\mu\nu}$"
	r0 = 0.5

	def plot(self, analysis_type = "bootstrap"):
		# Populating the plot values
		self.plot_values = []

		# Retrieving data deepending on analysis type we are choosing
		if analysis_type == "bootstrap":
			data = self.bootstrap_data
		elif analysis_type == "jackknife":
			data = self.jackknife_data
		elif analysis_type == "unanalyzed":
			data = self.unanalyzed_data
		else:
			raise KeyError("Analysis %s not recognized" % analysis_type)

		# Sorts data into a format specific for the plotting method
		for beta in sorted(data.keys()):
			values = {}
			values["beta"] = beta
			values["x"] = self.flow_time/self.r0**2*getLatticeSpacing(beta)**2
			values["y"] = -data[beta]["y"]*self.flow_time**2/64.0
			values["y_err"] = data[beta]["y_error"]*self.flow_time**2/64.0
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta],beta)
			values["color"] = self.colors[beta]
			self.plot_values.append(values)

		self._plot_core(fname_appendix = "_" + analysis_type)#,x_limits=x_limits,y_limits=y_limits)

	def plot_continiuum(self):
		# Retrieves t0 values
		self._get_t0(self.plot_values)

		# Fits to continiuum and retrieves values to be plotted
		x, y, y_std, ar0, t0, t0_err, a, b, a_err, b_err = self._linefit_to_continiuum()

		# Creates empty arrays for populating with data points from beta values and continiuum limit
		x_points = np.zeros(len(ar0)+1)
		y_points = np.zeros(len(ar0)+1)
		y_points_err = np.zeros(len(ar0)+1)

		# Populates arrays with first fitted element
		x_points[0] = x[0]
		y_points[0] = y[0]
		y_points_err[0] = y_std[0]

		# Populates arrays with beta data points
		x_points[1:] = ar0
		y_points[1:] = t0
		y_points_err[1:] = t0_err

		# Creates figure and plot window
		fig = plt.figure(self.dpi)
		ax = fig.add_subplot(111)

		# ax.axvline(0,linestyle="--",color="0",alpha=0.5)
		ax.errorbar(x_points[1:],y_points[1:],yerr=y_points_err[1:],fmt="o",color="0",ecolor="0")
		ax.errorbar(x_points[0],y_points[0],yerr=y_points_err[0],fmt="o",capthick=4,color="r",ecolor="r")
		ax.plot(x,y,color="0")#,label=r"$y=%2.4fx + %2.4f$" % (a,b))
		ax.set_ylabel(r"$\frac{\sqrt{8t_0}}{r_0}$")
		ax.set_xlabel(r"$(a/r_0)^2$")
		# ax.set_title(r"Continiuum limit reference scale: $t_{0,cont}=%2.4f\pm%g$" % ((self.r0*y_points[0])**2/8,(self.r0*y_points_err[0])**2/8))
		ax.set_xlim(-0.005,0.045)
		ax.set_ylim(0.92,0.98)
		ax.grid(True)

		# Fixes axis tick intervals
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.01))

		# Puts on some arrows at relevant points
		ax.annotate(r"$a=0.05$fm", xy=((0.01/self.r0)**2, end), xytext=((0.01/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")
		ax.annotate(r"$a=0.07$fm", xy=((0.07/self.r0)**2, end), xytext=((0.07/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")
		ax.annotate(r"$a=0.1$fm", xy=((0.1/self.r0)**2, end), xytext=((0.1/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")

		# ax.legend(loc="lower left")

		# Saves figure
		fname = os.path.join(self.output_folder_path,"post_analysis_%s_continiuum.png" % self.observable_name_compact)
		fig.savefig(fname,dpi=self.dpi)

		print "Figure created in %s" % fname
		plt.show()
		plt.close(fig)

	def _get_t0(self,data_values):
		self.t0_values = []

		# print "KEYS SHOULD BE BETA VALUES: ", self.bs_raw.keys()
		# exit(1)

		# Populates values to be plotted and 
		for values in data_values:
			t0_batch = {}
			t0_batch["beta"] = values["beta"]
			t0_batch["t0"],t0_batch["t0_err"] = fit_line(	values["x"],self.bs_raw[values["beta"]][self.observable_name_compact],0.3,0.015, 
															data_function_correction = lambda x : -x*self.flow_time**2/64.0,
															fit_function_modifier = lambda x : x*self.r0**2,
															plot_fit_window = True)
			# print "BETA %f:" % t0_batch["beta"], t0_batch["t0"],t0_batch["t0_err"]
			# t0_batch["t0"],t0_batch["t0_err"] = self._linefit_t0(values["x"],values["y"],values["y_err"],t0_batch["beta"])
			# print "BETA %f:" % t0_batch["beta"], t0_batch["t0"],t0_batch["t0_err"]
			t0_batch["a"] = getLatticeSpacing(t0_batch["beta"])

			# Adds to list of batch-values
			self.t0_values.append(t0_batch)

	def _linefit_to_continiuum(self):
		x_datapoints = np.asarray([val["a"] for val in self.t0_values])
		b = np.asarray([val["beta"] for val in self.t0_values])
		t0 = np.asarray([val["t0"] for val in self.t0_values])
		t0_err = np.asarray([val["t0_err"] for val in self.t0_values])

		# Reversing arrays, since increasing beta is decreasing lattice spacing and sets up y axis values, x is already the a-values
		x_datapoints = (x_datapoints[::-1] / self.r0)**2
		y_datapoints = np.sqrt(8*t0[::-1]) / self.r0
		y_datapoints_error = (8*t0_err[::-1]) / (np.sqrt(8*t0[::-1])) / self.r0

		# Fitting data
		pol, polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x_datapoints, y_datapoints,sigma=y_datapoints_error)
		# print pol,"\n",polcov
		# pol, polcov = np.polyfit(x_datapoints,y_datapoints,1,rcond=None,full=False,w=1.0/y_datapoints_error,cov=True)
		# print pol,"\n",polcov

		# Gets line properties
		a = pol[0]
		b = pol[1]
		a_err, b_err = np.sqrt(np.diag(polcov))

		# Sets up line fitted variables		
		x = np.linspace(0,x_datapoints[-1]*1.03,100)
		y = a * x + b
		y_std = a_err * x + b_err

		return x,y,y_std,x_datapoints,y_datapoints,y_datapoints_error, a, b, a_err, b_err

def main(args):
	"""
	Args should be post-analysis folder
	"""
	# Loads data from post analysis folder
	data = PostAnalysisDataReader(args[0])
	print "Retrieving data from folder: %s" % args[0]

	# Rewrites all of the data to a single file for sharing with giovanni
	data.write_batch_to_single_file()

	# # Plots topsus
	# topsus_analysis = TopSusPostAnalysis(data.data_observables["topsus"],data.data_batch_name,data.flow_time)
	# topsus_analysis = TopSusPostAnalysis(data,"topsus")
	# topsus_analysis.plot()

	# Retrofits the topsus for continiuum limit

	# Plots energy
	# energy_analysis = EnergyPostAnalysis(data.data_observables["energy"],data.flow_time,data.data_batch_name)
	energy_analysis = EnergyPostAnalysis(data,"energy")
	energy_analysis.plot()

	# Retrofits the energy for continiuum limit
	energy_analysis.plot_continiuum()

if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		args = sys.argv[1:]
	else:
		args = ["../output/post_analysis_data/data2"]
		# args = ["../output/post_analysis_data_old"]
	main(args)