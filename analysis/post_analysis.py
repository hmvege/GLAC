import numpy as np, matplotlib.pyplot as plt, sys, os, scipy.optimize as sciopt
from tools.folderreadingtools import check_folder
from statistics.line_fitting import fit_line_form_bootstrap, fit_line
from tools.postanalysisdatareader import PostAnalysisDataReader, getLatticeSpacing

# import tqdm
# from scipy.optimize import curve_fit

class PostAnalysis:
	"""
	Post analysis base class
	"""
	observable_name = "Observable"
	observable_name_compact = "obs"
	x_label = r""
	y_label = r""
	dpi=350
	size_labels = {	6.0  : r"$24^3 \times 48$",
					6.1  : r"$28^3 \times 56$",
					6.2  : r"$32^3 \times 64$",
					6.45 : r"$48^3 \times 96$"}
	r0 = 0.5

	def __init__(self,data,observable,base_output_folder="../figures/post_analysis"):
		# Retrieves relevant data values
		self.flow_time 			= data.flow_time
		self.unanalyzed_data 	= {beta:data.data_observables[observable][beta]["unanalyzed"] for beta in sorted(data.data_observables[observable].keys())} # Data should now be sorted by beta values
		self.bootstrap_data 	= {beta:data.data_observables[observable][beta]["bootstrap"] for beta in sorted(data.data_observables[observable].keys())}
		self.jackknife_data 	= {beta:data.data_observables[observable][beta]["jackknife"] for beta in sorted(data.data_observables[observable].keys())}
		self.bs_raw 			= data.bs_data_raw
		self.jk_raw 			= data.jk_data_raw
		self.ac_corrections		= data.ac_corrections

		self.NBoots = self.bs_raw[self.bs_raw.keys()[0]][self.observable_name_compact].shape[-1]

		# Small test to ensure that the number of bootstraps and number of different beta batches match
		err_msg = "Number of bootstraps do not match number of different beta values"
		assert sum([True for i in self.bs_raw.keys() if self.bs_raw[i][self.observable_name_compact].shape[-1] == self.NBoots]) == data.N_betas, err_msg

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

	def _check_plot_values(self):
		"""
		Checks if we have set the analysis data type yet.
		"""
		if not hasattr(self,"plot_values"):
			raise AttributeError("set_analysis_data_type() has not been set yet.")

	def set_analysis_data_type(self,analysis_data_type="bootstrap"):
		self.plot_values = []

		# Retrieving data depending on analysis type we are choosing
		if analysis_data_type == "bootstrap":
			data = self.bootstrap_data
		elif analysis_data_type == "jackknife":
			data = self.jackknife_data
		elif analysis_data_type == "unanalyzed":
			data = self.unanalyzed_data
		else:
			raise KeyError("Analysis %s not recognized" % analysis_data_type)

		# Makes it a global constant so it can be added in plot figure name
		self.analysis_data_type = analysis_data_type

		# Initiates plot values
		self._initiate_plot_values(data)

	def _initiate_plot_values(self,*args,**kwargs):
		print "Warning: default plot value initiater not completed."
		exit(1)

	def plot(self, x_limits = False, y_limits = False):
		"""
		Function for making a basic plot of all the different beta values together.
		"""
		print "Plotting %s for betas %s together" % (self.observable_name_compact,", ".join([str(b) for b in sorted(self.unanalyzed_data.keys())]))
		fig = plt.figure(dpi=self.dpi)
		ax = fig.add_subplot(111)

		self._check_plot_values()

		# Retrieves values to plot
		for value in self.plot_values:
			x = value["x"]
			y = value["y"]
			y_err = value["y_err"]
			ax.plot(x,y,"-",label=value["label"],color=value["color"])
			ax.fill_between(x, y-y_err, y+y_err,alpha=0.5, edgecolor='', facecolor=value["color"])

		# print self.flow_time[1:]**2*self._energy_continiuum(self.flow_time[1:])[0]
		# ax.plot(self.flow_time[1:]/self.r0**2,self.flow_time[1:]**2*self._energy_continiuum(self.flow_time[1:])[0],color="b")

		# Basic plotting commands
		ax.grid(True)
		ax.set_title(r"%s %s" % (self.observable_name, self.formula))
		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)
		ax.legend(loc="lower right")

		# Sets axes limits if provided
		if x_limits != False:
			ax.set_xlim(x_limits)
		if y_limits != False:
			ax.set_ylim(y_limits)

		# Saves and closes figure
		fname = os.path.join(self.output_folder_path,"post_analysis_%s_%s.png" % (self.observable_name_compact,self.analysis_data_type))
		plt.savefig(fname)
		print "Figure saved in %s" % fname
		# plt.show()
		plt.close(fig)

	def _get_beta_values_to_fit(self,fit_target, fit_interval, axis,
								fit_type = "bootstrap_fit",
								fit_function_modifier = lambda x : x,
								plot_fit_window = False):
		"""
		Retrieves a line fitted value at a target t0.
		Available fit_types:
			- bootstrap_fit (default)
			- linefit_data
			- nearest
		"""
		self.beta_fit = []

		self._check_plot_values()

		# Populates values to be plotted and 
		for values in self.plot_values:
			bfit = {}
			# Sets beta value for data
			bfit["beta"] = values["beta"]

			# Retrieves fit value as well as its error
			if fit_type == "bootstrap_fit":
				bfit["t0"], bfit["t0_err"] = fit_line_form_bootstrap(	values["x"],values["bs"],self.observable_name_compact,
																		values["beta"],fit_target,fit_interval, axis=axis,
																		fit_function_modifier = fit_function_modifier,
																		plot_fit_window = plot_fit_window)
			elif fit_type == "data_line_fit":
				bfit["t0"],bfit["t0_err"] = fit_line(	values["x"],values["y"],values["y_err"],
														self.observable_name_compact,values["beta"],
														fit_target,fit_interval, axis=axis,
														fit_function_modifier = fit_function_modifier,
														plot_fit_window = plot_fit_window)
			elif fit_type == "nearest_val_fit":
				raise NotImplementedError("'nearest_val_fit' not implemented as a fit type yet.")
			else:
				raise KeyError("No fit_type named %s. Options: 'bootstrap_fit', 'data_line_fit' or 'nearest_val_fit'" % fit_type)

			# Adds lattice spacing to fit
			bfit["a"] = getLatticeSpacing(bfit["beta"])

			# Adds to list of batch-values
			self.beta_fit.append(bfit)

	@staticmethod
	def _get_line_prop(poly, polycov, x_points, y_points):
		"""
		Small static function for curve_fit and polyfit methods.
		"""
		# Gets line properties
		a = poly[0]
		b = poly[1]
		a_err, b_err = np.sqrt(np.diag(polycov))

		# Sets up line fitted variables		
		x = np.linspace(0,x_points[-1]*1.03,1000)
		y = a * x + b
		y_std = a_err * x + b_err

		return x, y, y_std, a, b, a_err, b_err

	def _linefit_to_continiuum(self,x_points,y_points,y_points_error,fit_type="least_squares"):
		"""
		Fits a a set of values to continiuum.
		Args:
			x_points (numpy float array) : x points to data fit
			y_points (numpy float array) : y points to data fit
			y_points_error (numpy float array) : error of y points to data fit
			[optional] fit_type (str) : type of fit to perform. Options: 'curve_fit' (default), 'polyfit'
		"""
		# Fitting data
		if fit_type == "least_squares":
			pol, polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x_points, y_points,sigma=y_points_error)
			return self._get_line_prop(pol, polcov, x_points, y_points)
		elif fit_type == "polynomial":
			pol, polcov = np.polyfit(x_points,y_points,1,rcond=None,full=False,w=1.0/y_points_error,cov=True)
			return self._get_line_prop(pol,polcov)
		elif fit_type == "interpolation":
			return self.interpolate(x_points, y_points, y_points_error)
			return None
		else:
			raise KeyError("fit_type '%s' not recognized." % fit_type)


	# def _linefit_to_continiuum(self,fit_type="least_squares"):
	# 	"""
	# 	Fits values found in _get_beta_values_to_fit, and fits them to continium.
	# 	Args:
	# 		fit_type 	(str) : type of fit to be performed. Options: curve_fit (default), polyfit 
	# 	"""
	# 	# Retrieves values to fit
	# 	x_datapoints = np.asarray([val["a"] for val in self.beta_fit])[::-1]
	# 	t0 = np.asarray([val["t0"] for val in self.beta_fit])
	# 	t0_err = np.asarray([val["t0_err"] for val in self.beta_fit])

	# 	# Reversing arrays, since increasing beta is decreasing lattice spacing and sets up y axis values, x is already the a-values
	# 	x_datapoints = (x_datapoints[::-1] / self.r0)**2
	# 	y_datapoints = np.sqrt(8*t0[::-1]) / self.r0
	# 	y_datapoints_error = (8*t0_err[::-1]) / (np.sqrt(8*t0[::-1])) / self.r0

	# 	# Fitting data
	# 	if fit_type == "least_squares":
	# 		pol, polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x_datapoints, y_datapoints,sigma=y_datapoints_error)
	# 	elif fit_type == "polynomial":
	# 		pol, polcov = np.polyfit(x_datapoints,y_datapoints,1,rcond=None,full=False,w=1.0/y_datapoints_error,cov=True)
	# 	else:
	# 		raise KeyError("fit_type '%s' not recognized." % fit_type)

	# 	# Gets line properties
	# 	a = pol[0]
	# 	b = pol[1]
	# 	a_err, b_err = np.sqrt(np.diag(polcov))

	# 	# Sets up line fitted variables		
	# 	x = np.linspace(0,x_datapoints[-1]*1.03,100)
	# 	y = a * x + b
	# 	y_std = a_err * x + b_err

	# 	return x,y,y_std,x_datapoints,y_datapoints,y_datapoints_error, a, b, a_err, b_err

class TopSusPostAnalysis(PostAnalysis):
	observable_name = "Topological Susceptibility"
	observable_name_compact = "topsus"

	# Regular plot variables
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"
	formula = r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$"

	# Continiuum plot variables
	y_label_continiuum = r"$\chi^{1/4}[GeV]$"
	# x_label_continiuum = r"$a/{{r_0}^2}$"
	x_label_continiuum = r"$a[fm]$"

	def _initiate_plot_values(self,data):
		"""
		Function that sorts data into a format specific for the plotting method
		"""		
		for beta in sorted(data.keys()):
			values = {}
			values["beta"] = beta
			values["a"] = getLatticeSpacing(beta)
			values["x"] = getLatticeSpacing(beta)*np.sqrt(8*self.flow_time)
			values["y"] = data[beta]["y"]
			values["bs"] = self.bs_raw[beta][self.observable_name_compact]
			values["y_err"] = data[beta]["y_error"]
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta],beta)
			values["color"] = self.colors[beta]
			self.plot_values.append(values)

	def plot_continiuum(self,fit_target,fit_interval,fit_type):
		# Retrieves t0 values used to be used for continium fitting
		self._get_beta_values_to_fit(fit_target, fit_interval, axis = "x", fit_type = fit_type, plot_fit_window = False)

		# Builts plot variables
		a_lattice_spacings = np.asarray([val["a"] for val in self.beta_fit])[::-1]
		t_fit_points = np.asarray([val["t0"] for val in self.beta_fit])[::-1]
		t_fit_points_errors = np.asarray([val["t0_err"] for val in self.beta_fit])[::-1]

		# Initiates empty arrays for
		x_points = np.zeros(len(a_lattice_spacings)+1)
		y_points = np.zeros(len(a_lattice_spacings)+1)
		y_points_err = np.zeros(len(a_lattice_spacings)+1)

		# Populates with fit data
		x_points[1:] = a_lattice_spacings
		y_points[1:] = t_fit_points
		y_points_err[1:] = t_fit_points_errors

		# Fits to continiuum and retrieves values to be plotted
		x_line, y_line, y_line_std, a, b, a_err, b_err = self._linefit_to_continiuum(x_points[1:], y_points[1:], y_points_err[1:])

		# Populates arrays with first fitted element
		x_points[0] = x_line[0]
		y_points[0] = y_line[0]
		y_points_err[0] = y_line_std[0]

		# Creates figure and plot window
		fig = plt.figure(self.dpi)
		ax = fig.add_subplot(111)

		# ax.axvline(0,linestyle="--",color="0",alpha=0.5)
		ax.errorbar(x_points[1:],y_points[1:],yerr=y_points_err[1:],fmt="o",color="0",ecolor="0")
		ax.errorbar(x_points[0],y_points[0],yerr=y_points_err[0],fmt="o",capthick=4,color="r",ecolor="r",label=r"$\chi^{1/4}=%.2f\pm%.2f$" % (y_points[0],y_points_err[0]))
		ax.plot(x_line,y_line,color="0")#,label=r"$y=(%.3f\pm%.3f)x + %.4f\pm%.4f$" % (a,a_err,b,b_err))
		ax.fill_between(x_line, y_line-y_line_std, y_line+y_line_std,alpha=0.2, edgecolor='', facecolor="0")
		ax.set_ylabel(self.y_label_continiuum)
		ax.set_xlabel(self.x_label_continiuum)

		ax.set_title(r"Continiuum limit at: $\sqrt{8t_{flow}} = %.2f[fm]$" % (fit_target))

		ax.legend()
		ax.grid(True)

		# Saves figure
		fname = os.path.join(self.output_folder_path,"post_analysis_%s_continiuum%s_%s.png" % (self.observable_name_compact,str(fit_target).replace(".",""),fit_type.strip("_")))
		fig.savefig(fname,dpi=self.dpi)

		print "Continiuum plot of %s created in %s" % (self.observable_name.lower(),fname)
		# plt.show()
		plt.close(fig)

class EnergyPostAnalysis(PostAnalysis):
	observable_name = "Energy"
	observable_name_compact = "energy"

	# Regular plot variables
	y_label = r"$t^2\langle E\rangle$"
	x_label = r"$t/r_0^2$"
	formula = r"$\langle E\rangle = -\frac{1}{64V}F_{\mu\nu}^a{F^a}^{\mu\nu}$"

	# Continiuum plot variables
	x_label_continiuum = r"$(a/r_0)^2$"
	y_label_continiuum = r"$\frac{\sqrt{8t_0}}{r_0}$"

	def _energy_continiuum(self,t):
		"""
		Second order approximation of the energy.
		"""
		coupling = self._coupling_alpha(t)
		k1 = 1.0978
		mean_E = 3.0/(4.0*np.pi*t**2) * coupling * (1 + k1 * coupling)
		mean_E_error = 3.0/(4.0*np.pi*t**2)*coupling**3
		return mean_E, mean_E_error

	@staticmethod
	def _coupling_alpha(t):
		"""
		Running coupling constant.
		"""
		q = 1.0/np.sqrt(8*t)
		beta0 = 11.0
		LAMBDA = 0.34 # [GeV]
		alpha = 4*np.pi / (beta0 * np.log(q/LAMBDA**2))
		return alpha

	def _function_correction(self,x):
		"""
		Function that corrects the energy data
		"""
		return x*self.flow_time**2

	def _initiate_plot_values(self,data):
		# Sorts data into a format specific for the plotting method
		for beta in sorted(data.keys()):
			values = {}
			values["beta"] = beta
			values["a"] = getLatticeSpacing(beta)
			values["x"] = self.flow_time/self.r0**2*getLatticeSpacing(beta)**2
			values["y"] = self._function_correction(data[beta]["y"])
			values["bs"] = np.asarray([self._function_correction(self.bs_raw[beta][self.observable_name_compact][:,iBoot]) for iBoot in xrange(self.NBoots)]).T
			values["y_err"] = - self._function_correction(data[beta]["y_error"]) # negative since the minus sign will go away during linear error propagation
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta],beta)
			values["color"] = self.colors[beta]

			self.plot_values.append(values)

	def _linefit_to_continiuum(self,x_points,y_points,y_points_error,fit_type="least_squares"):
		"""
		Fits a a set of values to continiuum.
		Args:
			x_points (numpy float array) : x points to data fit
			y_points (numpy float array) : y points to data fit
			y_points_error (numpy float array) : error of y points to data fit
			[optional] fit_type (str) : type of fit to perform. Options: 'curve_fit' (default), 'polyfit'
		"""
		# Fitting data
		if fit_type == "least_squares":
			pol, polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x_points, y_points,sigma=y_points_error)
		elif fit_type == "polynomial":
			pol, polcov = np.polyfit(x_points,y_points,1,rcond=None,full=False,w=1.0/y_points_error,cov=True)
		else:
			raise KeyError("fit_type '%s' not recognized." % fit_type)

		# Gets line properties
		a = pol[0]
		b = pol[1]
		a_err, b_err = np.sqrt(np.diag(polcov))

		# Sets up line fitted variables		
		x = np.linspace(0,x_points[-1]*1.03,1000)
		y = a * x + b
		y_std = a_err * x + b_err

		return x, y, y_std, a, b, a_err, b_err

	def plot_continiuum(self,fit_target,fit_interval,fit_type,plot_arrows = [0.05,0.07,0.1], legend_location = "best"):
		# Retrieves t0 values used to be used for continium fitting
		self._get_beta_values_to_fit(fit_target,fit_interval,axis="y",
									fit_type = fit_type, 
									fit_function_modifier = lambda x : x*self.r0**2,
									plot_fit_window = False)

		a_lattice_spacings = np.asarray([val["a"] for val in self.beta_fit])[::-1]
		t_fit_points = np.asarray([val["t0"] for val in self.beta_fit])[::-1]
		t_fit_points_errors = np.asarray([val["t0_err"] for val in self.beta_fit])[::-1]

		# Initiates empty arrays for
		x_points = np.zeros(len(a_lattice_spacings)+1)
		y_points = np.zeros(len(a_lattice_spacings)+1)
		y_points_err = np.zeros(len(a_lattice_spacings)+1)

		# Populates with fit data
		x_points[1:] = (a_lattice_spacings / self.r0)**2 
		y_points[1:] = np.sqrt(8*(t_fit_points)) / self.r0
		y_points_err[1:] = (8*t_fit_points_errors) / (np.sqrt(8*t_fit_points)) / self.r0

		# Fits to continiuum and retrieves values to be plotted
		x_line, y_line, y_line_std, a, b, a_err, b_err = self._linefit_to_continiuum(x_points[1:], y_points[1:], y_points_err[1:])

		# Populates arrays with first fitted element
		x_points[0] = x_line[0]
		y_points[0] = y_line[0]
		y_points_err[0] = y_line_std[0]

		# Creates figure and plot window
		fig = plt.figure(self.dpi)
		ax = fig.add_subplot(111)

		# ax.axvline(0,linestyle="--",color="0",alpha=0.5)

		ax.errorbar(x_points[1:],y_points[1:],yerr=y_points_err[1:],fmt="o",color="0",ecolor="0")
		ax.errorbar(x_points[0],y_points[0],yerr=y_points_err[0],fmt="o",capthick=4,color="r",ecolor="r")
		ax.plot(x_line,y_line,color="0",label=r"$y=(%.3f\pm%.3f)x + %.4f\pm%.4f$" % (a,a_err,b,b_err))
		ax.fill_between(x_line, y_line-y_line_std, y_line+y_line_std,alpha=0.2, edgecolor='', facecolor="0")
		ax.set_ylabel(self.y_label_continiuum)
		ax.set_xlabel(self.x_label_continiuum)

		# ax.set_title(r"Continiuum limit reference scale: $t_{0,cont}=%2.4f\pm%g$" % ((self.r0*y_points[0])**2/8,(self.r0*y_points_err[0])**2/8))

		ax.set_xlim(-0.005,0.045)
		ax.set_ylim(0.92,0.98)
		ax.legend()
		ax.grid(True)

		# Fixes axis tick intervals
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.01))

		# Puts on some arrows at relevant points
		for arrow in plot_arrows:
			ax.annotate(r"$a=%.2g$fm" % arrow, xy=((arrow/self.r0)**2, end), xytext=((arrow/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")
		
		ax.legend(loc=legend_location) # "lower left"

		# Saves figure
		fname = os.path.join(self.output_folder_path,"post_analysis_%s_continiuum%s_%s.png" % (self.observable_name_compact,str(fit_target).replace(".",""),fit_type.strip("_")))
		fig.savefig(fname,dpi=self.dpi)

		print "Continiuum plot of %s created in %s" % (self.observable_name.lower(),fname)
		plt.close(fig)

	def coupling_fit(self):
		print "Finding Lambda"


		pass

def main(args):
	"""
	Args should be post-analysis folder
	"""
	# Loads data from post analysis folder
	data = PostAnalysisDataReader(args[0])
	print "Retrieving data from folder: %s" % args[0]

	# Rewrites all of the data to a single file for sharing with giovanni
	# data.write_batch_to_single_file()

	# Plots topsus
	topsus_analysis = TopSusPostAnalysis(data,"topsus")
	topsus_analysis.set_analysis_data_type("bootstrap")
	topsus_analysis.plot()

	# Retrofits the topsus for continuum limit
	continium_targets = [0.3,0.4,0.5,0.58]
	for cont_target in continium_targets:
		topsus_analysis.plot_continiuum(cont_target,0.015,"data_line_fit")

	# Plots energy
	energy_analysis = EnergyPostAnalysis(data,"energy")
	energy_analysis.set_analysis_data_type("bootstrap")
	energy_analysis.plot()

	# Retrofits the energy for continiuum limit
	energy_analysis.plot_continiuum(0.3, 0.015,"bootstrap_fit")

	# Plot running coupling
	energy_analysis.coupling_fit()

if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		args = sys.argv[1:]
	else:
		args = ["../output/post_analysis_data/data2"]
		# args = ["../output/post_analysis_data/data4"]
	main(args)