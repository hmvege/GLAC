import numpy as np, matplotlib.pyplot as plt, sys, os, scipy.optimize as sciopt
from tools.folderreadingtools import check_folder
# from statistics.line_fitting import fit_line_form_bootstrap, fit_line
from statistics.linefit import LineFit
from tools.postanalysisdatareader import PostAnalysisDataReader, getLatticeSpacing
import itertools

# import tqdm
# from scipy.optimize import curve_fit

__all__ = [
	"EnergyPostAnalysis", "PlaqPostAnalysis", 
	# Different topological susceptibility definitions
	"TopSusPostAnalysis", "QtQ0PostAnalysis", "TopCharge4PostAnalysis",
	"TopSusTPostAnalysis", "TopSusEuclSplitPostAnalysis",
	"TopSusMCSplitPostAnalysis",
	# Different topological charge definitions
	"TopChargePostAnalysis", "TopChargeTPostAnalysis", 
	"TopChargeEuclSplitPostAnalysis", "TopChargeMCSplitPostAnalysis",
]

class _PostAnalysis:
	"""Post analysis base class."""
	observable_name = "Observable"
	observable_name_compact = "obs"
	formula = ""
	x_label = r""
	y_label = r""
	dpi=350
	size_labels = {	6.0  : r"$24^3 \times 48$",
					6.1  : r"$28^3 \times 56$",
					6.2  : r"$32^3 \times 64$",
					6.45 : r"$48^3 \times 96$"}
	r0 = 0.5
	sub_obs = False
	observable_intervals = {}
	# blue, green, red purple
	beta_colors = ["#5cbde0", "#6fb718", "#bc232e", "#8519b7"]

	def __init__(self, data, with_autocorr=True, figures_folder="../figures", verbose=False):
		if with_autocorr:
			self.ac = "with_autocorr"
		else:
			self.ac = "without_autocorr"
		self.with_autocorr = with_autocorr
		observable = self.observable_name_compact

		self.verbose = verbose

		# Retrieves relevant data values and sorts them by beta values
		self.flow_time = data.flow_time
		self.unanalyzed_data = {}
		self.bootstrap_data	= {}
		self.jackknife_data = {}

		for beta in sorted(data.data_observables[observable].keys()):
			if self.sub_obs:
				self.observable_intervals[beta] = data.data_observables[observable][beta].keys()
				if not beta in self.unanalyzed_data:
					self.unanalyzed_data[beta] = {}
					self.bootstrap_data[beta] = {}
					self.jackknife_data[beta] = {}
				for subobs in data.data_observables[observable][beta]:
					self.unanalyzed_data[beta][subobs] = data.data_observables[observable][beta][subobs][self.ac]["unanalyzed"]
					self.bootstrap_data[beta][subobs] = data.data_observables[observable][beta][subobs][self.ac]["bootstrap"]
					self.jackknife_data[beta][subobs] = data.data_observables[observable][beta][subobs][self.ac]["jackknife"]
			else:
				self.unanalyzed_data[beta] = data.data_observables[observable][beta][self.ac]["unanalyzed"]
				self.bootstrap_data[beta] = data.data_observables[observable][beta][self.ac]["bootstrap"]
				self.jackknife_data[beta] = data.data_observables[observable][beta][self.ac]["jackknife"]

		self.bs_raw = data.raw_analysis["bootstrap"]
		self.jk_raw = data.raw_analysis["jackknife"]
		self.ac_corrections	= data.raw_analysis["autocorrelation"]

		if not isinstance(self.bs_raw[self.bs_raw.keys()[0]][self.observable_name_compact], np.ndarray):
			self.NBoots = self.bs_raw.values()[0][self.observable_name_compact].values()[0].shape[-1]
		else:
			self.NBoots = self.bs_raw[self.bs_raw.keys()[0]][self.observable_name_compact].shape[-1]

		# Small test to ensure that the number of bootstraps and number of different beta batches match
		err_msg = "Number of bootstraps do not match number of different beta values"
		if self.sub_obs:
			chk_bs_len = lambda _a, _i: _a[_i][self.observable_name_compact].values()[0].shape[-1] == self.NBoots
		else:
			chk_bs_len = lambda _a, _i: _a[_i][self.observable_name_compact].shape[-1] == self.NBoots

		assert sum([True for i in self.bs_raw.keys() if chk_bs_len(self.bs_raw, i)]) == data.N_betas, err_msg

		# Creates base output folder for post analysis figures
		self.figures_folder = figures_folder
		check_folder(self.figures_folder, dryrun=False, verbose=self.verbose)
		check_folder(os.path.join(self.figures_folder, data.batch_name), dryrun=False, verbose=self.verbose)

		# Creates output folder
		self.post_anlaysis_folder = os.path.join(self.figures_folder, data.batch_name, "post_analysis")
		check_folder(self.post_anlaysis_folder, dryrun=False, verbose=self.verbose)

		# Creates observable output folder
		self.output_folder_path = os.path.join(self.post_anlaysis_folder, self.observable_name_compact)
		check_folder(self.output_folder_path, dryrun=False, verbose=self.verbose)

		# Creates colors to use
		self.colors = {}
		for color, beta in zip(self.beta_colors, sorted(data.data_observables[observable].keys())):
			self.colors[beta] = color

	def _check_plot_values(self):
		"""Checks if we have set the analysis data type yet."""
		if not hasattr(self,"plot_values"):
			raise AttributeError("set_analysis_data_type() has not been set yet.")

	def _get_analysis_data(self, analysis_data_type):
		"""Retrieving data depending on analysis type we are choosing"""
		if analysis_data_type == "bootstrap":
			return self.bootstrap_data
		elif analysis_data_type == "jackknife":
			return self.jackknife_data
		elif analysis_data_type == "unanalyzed":
			return self.unanalyzed_data
		else:
			raise KeyError("Analysis %s not recognized" % analysis_data_type)

	def set_analysis_data_type(self, analysis_data_type="bootstrap"):
		self.plot_values = {}

		data = self._get_analysis_data(analysis_data_type)

		# Makes it a global constant so it can be added in plot figure name
		self.analysis_data_type = analysis_data_type

		# Initiates plot values
		self._initiate_plot_values(data)

	def _initiate_plot_values(self, data):
		# Sorts data into a format specific for the plotting method
		for beta in sorted(data.keys()):
			if beta == 6.45: self.flow_time *= 2
			values = {}
			values["a"] = getLatticeSpacing(beta)
			values["x"] = values["a"]* np.sqrt(8*self.flow_time)
			values["y"] = data[beta]["y"]
			values["bs"] = self.bs_raw[beta][self.observable_name_compact]
			values["y_err"] = data[beta]["y_error"] # negative since the minus sign will go away during linear error propagation
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta], beta)
			values["color"] = self.colors[beta]
			self.plot_values[beta] = values

	def plot(self, x_limits=False, y_limits=False, plot_with_formula=False):
		"""
		Function for making a basic plot of all the different beta values together.

		Args:
			x_limits: limits of the x-axis. Default is False.
			y_limits: limits of the y-axis. Default is False.
			plot_with_formula: bool, default is false, is True will look for 
				formula for the y-value to plot in title.
		"""
		if self.verbose:
			print "Plotting %s for betas %s together" % (
				self.observable_name_compact,
				", ".join([str(b) for b in sorted(self.unanalyzed_data.keys())]))

		fig = plt.figure(dpi=self.dpi)
		ax = fig.add_subplot(111)

		self._check_plot_values()

		# Retrieves values to plot
		for beta in sorted(self.plot_values):
			value = self.plot_values[beta]
			x = value["x"]
			y = value["y"]
			y_err = value["y_err"]
			ax.plot(x, y, "-", label=value["label"], color=value["color"])
			ax.fill_between(x, y - y_err, y + y_err, alpha=0.5, edgecolor='', facecolor=value["color"])

		# print self.flow_time[1:]**2*self._energy_continuum(self.flow_time[1:])[0]
		# ax.plot(self.flow_time[1:]/self.r0**2,self.flow_time[1:]**2*self._energy_continuum(self.flow_time[1:])[0],color="b")

		# Sets the title string
		title_string = r"%s" % self.observable_name
		if plot_with_formula:
			title_string += r" %s" % self.formula

		# Basic plotting commands
		ax.grid(True)
		ax.set_title(title_string)
		ax.set_xlabel(self.x_label)
		ax.set_ylabel(self.y_label)
		ax.legend(loc="best", prop={"size": 8})

		# Sets axes limits if provided
		if x_limits != False:
			ax.set_xlim(x_limits)
		if y_limits != False:
			ax.set_ylim(y_limits)

		plt.tight_layout()

		# Saves and closes figure
		fname = self._get_plot_figure_name()
		plt.savefig(fname)
		print "Figure saved in %s" % fname
		# plt.show()
		plt.close(fig)

	def _get_beta_values_to_fit(self, fit_target, fit_interval, axis,
								fit_type="bootstrap_fit",
								fit_function_modifier=lambda x: x,
								plot_fit_window=False):
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
		for beta in self.plot_values:
			bfit = {}
			# Sets beta value for data
			bfit["beta"] = beta

			# Retrieves fit value as well as its error
			if fit_type == "bootstrap_fit":
				bfit["t0"], bfit["t0_err"] = fit_line_from_bootstrap(	values["x"], values["bs"], self.observable_name_compact,
																		beta, fit_target, fit_interval, axis=axis,
																		fit_function_modifier=fit_function_modifier,
																		plot_fit_window=plot_fit_window)
			elif fit_type == "data_line_fit":
				bfit["t0"], bfit["t0_err"] = fit_line(	values["x"], values["y"], values["y_err"],
														self.observable_name_compact, beta,
														fit_target, fit_interval, axis=axis,
														fit_function_modifier=fit_function_modifier,
														plot_fit_window=plot_fit_window)
			elif fit_type == "nearest_val_fit":
				raise NotImplementedError("'nearest_val_fit' not implemented as a fit type yet.")
			else:
				raise KeyError("No fit_type named %s. Options: 'bootstrap_fit', 'data_line_fit' or 'nearest_val_fit'" % fit_type)

			# Adds lattice spacing to fit
			bfit["a"] = getLatticeSpacing(bfit["beta"])

			# Adds to list of batch-values
			self.beta_fit.append(bfit)

	def _get_plot_figure_name(self):
		"""Retrieves appropriate figure file name."""
		output_folder = self.output_folder_path
		fname = "post_analysis_%s_%s.png" % (self.observable_name_compact, self.analysis_data_type)
		return os.path.join(output_folder, fname)

	def __str__(self):
		msg = "\n" +"="*100
		msg += "\nPost analaysis for:        " + self.observable_name_compact
		msg += "\nIncluding autocorrelation: " + self.ac
		msg += "\nOutput folder:             " + self.output_folder_path
		msg += "\n" + "="*100
		return msg

class _MultiObservablePostAnalysis(_PostAnalysis):
	"""
	Class to be inheritedfrom in case we got intervals or sub elements of the 
	same observable.
	"""
	sub_obs = True
	analysis_data_type = "bootstrap"

	def _initiate_plot_values(self, data, interval_index=None):
		# Sorts data into a format specific for the plotting method
		for beta in sorted(data.keys()):
			values = {}
			if interval_index == None:
				# Case where we have sub sections of observables, e.g. in euclidean time
				for sub_obs in self.observable_intervals[beta]:
					if beta == 6.45: self.flow_time *= 2
					sub_values = {}
					sub_values["a"] = getLatticeSpacing(beta)
					sub_values["x"] = sub_values["a"]* np.sqrt(8*self.flow_time)
					sub_values["y"] = data[beta][sub_obs]["y"]
					sub_values["bs"] = self.bs_raw[beta][self.observable_name_compact][sub_obs]
					sub_values["y_err"] = data[beta][sub_obs]["y_error"] # negative since the minus sign will go away during linear error propagation
					sub_values["label"] = r"%s $\beta=%2.2f$ %s" % (self.size_labels[beta], beta, sub_obs)
					sub_values["color"] = self.colors[beta]
					values[sub_obs] = sub_values
			else:
				sorted_intervals = sorted(data[beta].keys())
				values["a"] = getLatticeSpacing(beta)
				values["x"] = values["a"]* np.sqrt(8*self.flow_time)
				values["y"] = data[beta][sorted_intervals[interval_index]]["y"]
				values["bs"] = self.bs_raw[beta][self.observable_name_compact][sorted_intervals[interval_index]]
				values["y_err"] = data[beta][sorted_intervals[interval_index]]["y_error"] # negative since the minus sign will go away during linear error propagation
				values["label"] = r"%s $\beta=%2.2f$ %s" % (self.size_labels[beta], beta, sorted_intervals[interval_index])
				values["color"] = self.colors[beta]
				values["interval"] = sorted_intervals[interval_index]
			self.plot_values[beta] = values

	def set_analysis_type(self, analysis_data_type):
		"""Sets a global analysis type."""
		self.analysis_data_type = analysis_data_type

	def plot_interval(self, interval_index, **kwargs):
		"""Sets and plots only one interval."""
		self.interval_index = interval_index
		self.plot_values = {}
		data = self._get_analysis_data(self.analysis_data_type)
		self._initiate_plot_values(data, interval_index=interval_index)
		# Makes it a global constant so it can be added in plot figure name
		self.plot(**kwargs)

	def _get_plot_figure_name(self):
		"""Retrieves appropriate figure file name."""
		output_folder = os.path.join(self.output_folder_path, "slices")
		check_folder(output_folder, False, True)
		fname = "post_analysis_%s_%s_int%d.png" % (self.observable_name_compact, self.analysis_data_type, self.interval_index)
		return os.path.join(output_folder, fname)

	def get_N_intervals(self):
		"""Returns possible intervals for us to plot."""
		if self.verbose:
			print "Intervals N=%d, possible for %s: " % (len(self.observable_intervals), 
				self.observable_name_compact),
			print self.observable_intervals
		return len(self.observable_intervals), self.observable_intervals

	def plot_series(self, indexes, beta="all", x_limits=False, 
		y_limits=False, plot_with_formula=False):
		"""
		Method for plotting 4 axes together.

		Args:
			indexes: list containing integers of which intervals to plot together.
			beta: beta values to plot. Default is "all". Otherwise, 
				a list of numbers or a single beta value is provided.
			x_limits: limits of the x-axis. Default is False.
			y_limits: limits of the y-axis. Default is False.
			plot_with_formula: bool, default is false, is True will look for 
				formula for the y-value to plot in title.
		"""
		self.plot_values = {}
		data = self._get_analysis_data(self.analysis_data_type)
		self._initiate_plot_values(data)

		old_rc_paramx = plt.rcParams['xtick.labelsize']
		old_rc_paramy = plt.rcParams['ytick.labelsize']
		plt.rcParams['xtick.labelsize'] = 6
		plt.rcParams['ytick.labelsize'] = 6

		# Starts plotting
		# fig = plt.figure(sharex=True)
		fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)

		# Ensures beta is a list
		if not isinstance(beta, list):
			beta = [beta]

		# Sets the beta values to plot
		if beta[0] == "all" and len(beta) == 1:
			bvalues = self.plot_values
		else:
			bvalues = beta

		# print axes
		for ax, i in zip(list(itertools.chain(*axes)), indexes):
			for ibeta in bvalues:
				# Retrieves the values deepending on the indexes provided and beta values
				value = self.plot_values[ibeta][sorted(self.observable_intervals[ibeta])[i]]
				x = value["x"]
				y = value["y"]
				y_err = value["y_err"]
				ax.plot(x, y, "-", label=value["label"], color=value["color"])
				ax.fill_between(x, y - y_err, y + y_err, alpha=0.5, edgecolor='', facecolor=value["color"])
				
				# Basic plotting commands
				ax.grid(True)
				ax.legend(loc="best", prop={"size":5})

				# Sets axes limits if provided
				if x_limits != False:
					ax.set_xlim(x_limits)
				if y_limits != False:
					ax.set_ylim(y_limits)

		# Set common labels
		# https://stackoverflow.com/questions/6963035/pyplot-axes-labels-for-subplots
		fig.text(0.52, 0.035, self.x_label, ha='center', va='center', fontsize=9)
		fig.text(0.03, 0.5, self.y_label, ha='center', va='center', rotation='vertical', fontsize=11)

		# Sets the title string
		title_string = r"%s" % self.observable_name
		if plot_with_formula:
			title_string += r" %s" % self.formula
		plt.suptitle(title_string)
		plt.tight_layout(pad=1.7)

		# Saves and closes figure
		if beta == "all":
			folder_name = "beta%s" % beta
		else:
			folder_name = "beta%s" % "-".join([str(i) for i in beta])
		folder_name += "_N%s" % "".join([str(i) for i in indexes])
		folder_path = os.path.join(self.output_folder_path, folder_name)
		check_folder(folder_path, False, True)

		fname = os.path.join(folder_path, "post_analysis_%s_%s.png" % (self.observable_name_compact, self.analysis_data_type))
		plt.savefig(fname, dpi=400)
		print "Figure saved in %s" % fname
		# plt.show()
		plt.close(fig)

		plt.rcParams['xtick.labelsize'] = old_rc_paramx
		plt.rcParams['ytick.labelsize'] = old_rc_paramy


class TopSusPostAnalysis(_PostAnalysis):
	observable_name = "Topological Susceptibility"
	observable_name_compact = "topsus"

	# Regular plot variables
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi_t^{1/4}[GeV]$"
	formula = r"$\chi_t^{1/4}=\frac{\hbar c}{aV^{1/4}}\langle Q^2 \rangle^{1/4}$"

	# Continuum plot variables
	y_label_continuum = r"$\chi^{1/4}[GeV]$"
	# x_label_continuum = r"$a/{{r_0}^2}$"
	x_label_continuum = r"$a[fm]$"

	def _initiate_plot_values(self, data):
		"""
		Function that sorts data into a format specific for the plotting method
		"""		
		for beta in sorted(data.keys()):
			if beta == 6.45: self.flow_time *= 2
			values = {}
			values["a"] = getLatticeSpacing(beta)
			values["x"] = values["a"]*np.sqrt(8*self.flow_time)
			values["y"] = data[beta]["y"]
			values["bs"] = self.bs_raw[beta][self.observable_name_compact]
			values["y_err"] = data[beta]["y_error"]
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta], beta)
			values["color"] = self.colors[beta]
			self.plot_values[beta] = values

	def plot_continuum(self, fit_target):
		# Gets the beta values

		a, obs, obs_err = [], [], []
		for beta in sorted(self.plot_values):
			x = self.plot_values[beta]["x"]
			y = self.plot_values[beta]["y"]
			y_err = self.plot_values[beta]["y_err"]

			fit_index = np.argmin(np.abs(x - fit_target))

			a.append(self.plot_values[beta]["a"])
			obs.append(y[fit_index])
			obs_err.append(y_err[fit_index])

		# Initiates empty arrays for the continuum limit
		a = np.asarray(a)[::-1]
		obs = np.asarray(obs)[::-1]
		obs_err = np.asarray(obs_err)[::-1]

		# Continuum limit arrays
		N_cont = 1000
		a_cont = np.linspace(-0.01, a[-1]*1.1, N_cont)

		# Fits to continuum and retrieves values to be plotted
		continuum_fit = LineFit(a, obs, obs_err)

		y_cont, y_cont_err, fit_params, chi_squared = continuum_fit.fit_weighted(a_cont)

		# continuum_fit.plot(True)
		
		# Gets the continium value and its error
		cont_index = np.argmin(np.abs(a_cont))
		a0 = [a_cont[cont_index], a_cont[cont_index]]
		y0 = [y_cont[cont_index], y_cont[cont_index]]

		if y_cont_err[1][cont_index] < y_cont_err[0][cont_index]:
			y0_err_lower = y0[0] - y_cont_err[1][cont_index]
			y0_err_upper = y_cont_err[0][cont_index] - y0[0]
		else:
			y0_err_lower = y0[0] - y_cont_err[0][cont_index]
			y0_err_upper = y_cont_err[1][cont_index] - y0[0]

		y0_err = [[y0_err_lower, 0], [y0_err_upper, 0]]

		# Creates figure and plot window
		fig = plt.figure()
		ax = fig.add_subplot(111)

		# Plots linefit with errorband
		ax.plot(a_cont, y_cont, color="tab:blue", alpha=0.5)
		ax.fill_between(a_cont, y_cont_err[0], y_cont_err[1],
			alpha=0.5, edgecolor='', facecolor="tab:blue")

		# Plot lattice points
		ax.errorbar(a, obs, yerr=obs_err, fmt="o",
			color="tab:orange", ecolor="tab:orange")

		# plots continuum limit
		ax.errorbar(a0, y0, yerr=y0_err, fmt="o", capsize=None, # 5 is a good value for cap size
			capthick=1, color="tab:red", ecolor="tab:red",
			label=r"$\chi^{1/4}=%.3f\pm%.3f$" % (y0[0], (y0_err_lower + y0_err_upper)/2.0))

		ax.set_ylabel(self.y_label_continuum)
		ax.set_xlabel(self.x_label_continuum)
		ax.set_title(r"$\sqrt{8t_{flow,0}} = %.2f[fm], \chi^2 = %.2g$" % (fit_target, chi_squared))
		ax.set_xlim(-0.01, a[-1]*1.1)
		ax.legend()
		ax.grid(True)

		# print "Target: %.16f +/- %.16f" % (c[0], y0_err[0])

		# Saves figure
		fname = os.path.join(self.output_folder_path, 
			"post_analysis_%s_continuum%s_%s.png" % (self.observable_name_compact, 
				str(fit_target).replace(".",""), self.analysis_data_type))
		fig.savefig(fname, dpi=self.dpi)

		print "Continuum plot of %s created in %s" % (self.observable_name.lower(), fname)
		# plt.show()
		plt.close(fig)


class EnergyPostAnalysis(_PostAnalysis):
	observable_name = "Energy"
	observable_name_compact = "energy"

	# Regular plot variables
	y_label = r"$t^2\langle E\rangle$"
	x_label = r"$t/r_0^2$"
	formula = r"$\langle E\rangle = -\frac{1}{64V}F_{\mu\nu}^a{F^a}^{\mu\nu}$"

	# Continuum plot variables
	x_label_continuum = r"$(a/r_0)^2$"
	y_label_continuum = r"$\frac{\sqrt{8t_0}}{r_0}$"

	def _get_fit_interval(self, fit_target, fit_interval, y_mean):
		"""Function for finding the fit interval."""
		start_index = np.argmin(np.abs(y_mean - (fit_target - fit_interval)))
		end_index = np.argmin(np.abs(y_mean - (fit_target + fit_interval)))
		return start_index, end_index

	def _inverse_beta_fit(self, fit_target, fit_interval):
		"""
		Perform an inverse fit on the observable susceptibility and extracts 
		extracting x0 with x0 error.
		"""
		for beta in sorted(self.plot_values):
			x = self.plot_values[beta]["x"]
			y = self.plot_values[beta]["y"]
			y_err = self.plot_values[beta]["y_err"]

			index_low, index_high = self._get_fit_interval(fit_target, fit_interval, x)

			x = x[index_low:index_high]
			y = y[index_low:index_high]
			y_err = y_err[index_low:index_high]

			fit = LineFit(x, y, y_err)
			y_hat, y_hat_err, fit_params, chi_squared = fit.fit_weighted()
			b0, b0_err, b1, b1_err = fit_params

			self.plot_values[beta]["fit"] = {
				"y_hat": y_hat,
				"y_hat_err": y_hat_err,
				"b0": b0,
				"b0_err": b0_err,
				"b1": b1,
				"b1_err": b1_err,
				"chi_squared": chi_squared,
			}

			x0, x0_err = fit.inverse_fit(fit_target, weigthed=True)

			fit.plot(True)

			self.plot_values[beta]["fit"]["inverse"] = {
					"x0": x0,
					"x0_err": x0_err, 
			}

	def _energy_continuum(self, t):
		"""
		Second order approximation of the energy.
		"""
		coupling = self._coupling_alpha(t)
		k1 = 1.0978
		mean_E = 3.0/(4.0*np.pi*t**2) * coupling * (1 + k1*coupling)
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

	def _function_correction(self, x):
		"""Function that corrects the energy data."""
		return x*self.flow_time**2

	def _initiate_plot_values(self, data):
		# Sorts data into a format specific for the plotting method
		for beta in sorted(data.keys()):
			values = {}
			values["beta"] = beta
			values["a"] = getLatticeSpacing(beta)
			values["x"] = self.flow_time/self.r0**2*getLatticeSpacing(beta)**2
			values["y"] = self._function_correction(data[beta]["y"])
			values["bs"] = np.asarray([self._function_correction(self.bs_raw[beta][self.observable_name_compact][:,iBoot]) for iBoot in xrange(self.NBoots)]).T
			values["y_err"] = - self._function_correction(data[beta]["y_error"]) # negative since the minus sign will go away during linear error propagation
			values["label"] = r"%s $\beta=%2.2f$" % (self.size_labels[beta], beta)
			values["color"] = self.colors[beta]

			self.plot_values[beta] = values

	def _linefit_to_continuum(self, x_points, y_points, y_points_error, fit_type="least_squares"):
		"""
		Fits a a set of values to continuum.
		Args:
			x_points (numpy float array) : x points to data fit
			y_points (numpy float array) : y points to data fit
			y_points_error (numpy float array) : error of y points to data fit
			[optional] fit_type (str) : type of fit to perform. Options: 'curve_fit' (default), 'polyfit'
		"""

		# Fitting data
		if fit_type == "least_squares":
			pol, polcov = sciopt.curve_fit(lambda x, a, b: x*a + b, x_points, y_points, sigma=y_points_error)
		elif fit_type == "polynomial":
			pol, polcov = np.polyfit(x_points, y_points, 1, rcond=None, full=False, w=1.0/y_points_error, cov=True)
		else:
			raise KeyError("fit_type '%s' not recognized." % fit_type)

		# Gets line properties
		a = pol[0]
		b = pol[1]
		a_err, b_err = np.sqrt(np.diag(polcov))

		# Sets up line fitted variables		
		x = np.linspace(0, x_points[-1]*1.03, 1000)
		y = a * x + b
		y_std = a_err*x + b_err

		return x, y, y_std, a, b, a_err, b_err

	def plot_continuum(self, fit_target, fit_interval, fit_type, plot_arrows=[0.05, 0.07, 0.1], legend_location="best", line_fit_type="least_squares"):
		# Retrieves t0 values used to be used for continium fitting
		self._get_beta_values_to_fit(fit_target, fit_interval, axis="y",
									fit_type=fit_type, 
									fit_function_modifier=lambda x: x*self.r0**2,
									plot_fit_window=False)

		a_lattice_spacings = np.asarray([val["a"] for val in self.beta_fit])[::-1]
		t_fit_points = np.asarray([val["t0"] for val in self.beta_fit])[::-1]
		t_fit_points_errors = np.asarray([val["t0_err"] for val in self.beta_fit])[::-1]

		# Initiates empty arrays for
		x_points = np.zeros(len(a_lattice_spacings) + 1)
		y_points = np.zeros(len(a_lattice_spacings) + 1)
		y_points_err = np.zeros(len(a_lattice_spacings) + 1)

		# Populates with fit data
		x_points[1:] = (a_lattice_spacings / self.r0)**2 
		y_points[1:] = np.sqrt(8*(t_fit_points)) / self.r0
		y_points_err[1:] = (8*t_fit_points_errors) / (np.sqrt(8*t_fit_points)) / self.r0

		# Fits to continuum and retrieves values to be plotted
		x_line, y_line, y_line_std, a, b, a_err, b_err = self._linefit_to_continuum(x_points[1:], y_points[1:], y_points_err[1:], fit_type=line_fit_type)

		# Populates arrays with first fitted element
		x_points[0] = x_line[0]
		y_points[0] = y_line[0]
		y_points_err[0] = y_line_std[0]

		# Creates figure and plot window
		fig = plt.figure(self.dpi)
		ax = fig.add_subplot(111)

		# ax.axvline(0,linestyle="--",color="0",alpha=0.5)

		ax.errorbar(x_points[1:], y_points[1:], yerr=y_points_err[1:], fmt="o", color="0", ecolor="0")
		ax.errorbar(x_points[0], y_points[0], yerr=y_points_err[0], fmt="o", capthick=4, color="r", ecolor="r")
		ax.plot(x_line, y_line, color="0", label=r"$y=(%.3f\pm%.3f)x + %.4f\pm%.4f$" % (a, a_err, b, b_err))
		ax.fill_between(x_line, y_line-y_line_std, y_line + y_line_std, alpha=0.2, edgecolor='', facecolor="0")
		ax.set_ylabel(self.y_label_continuum)
		ax.set_xlabel(self.x_label_continuum)

		# ax.set_title(r"Continuum limit reference scale: $t_{0,cont}=%2.4f\pm%g$" % ((self.r0*y_points[0])**2/8,(self.r0*y_points_err[0])**2/8))

		ax.set_xlim(-0.005, 0.045)
		ax.set_ylim(0.92, 0.98)
		ax.legend()
		ax.grid(True)

		# Fixes axis tick intervals
		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.01))

		# Puts on some arrows at relevant points
		for arrow in plot_arrows:
			ax.annotate(r"$a=%.2g$fm" % arrow, xy=((arrow/self.r0)**2, end), xytext=((arrow/self.r0)**2, end+0.005), arrowprops=dict(arrowstyle="->"), ha="center")
		
		ax.legend(loc=legend_location) # "lower left"

		# Saves figure
		fname = os.path.join(self.output_folder_path, "post_analysis_%s_continuum%s_%s.png" % (self.observable_name_compact, str(fit_target).replace(".",""), fit_type.strip("_")))
		fig.savefig(fname, dpi=self.dpi)

		print "Continuum plot of %s created in %s" % (self.observable_name.lower(), fname)
		plt.close(fig)

	def coupling_fit(self):
		print "Finding Lambda"

		pass

class PlaqPostAnalysis(_PostAnalysis):
	"""Post-analysis of the topological charge."""
	observable_name = "Plaquette"
	observable_name_compact = "plaq"
	y_label = r"$P$"
	x_label = r"$\sqrt{8t}$[fm]"
	formula = r"$P = \frac{1}{16V} \sum_{x,\mu,\nu} \tr\mathcal{Re} P_{\mu\nu}$"

#### Topological charge definitions ####
class TopChargePostAnalysis(_PostAnalysis):
	"""Post-analysis of the topological charge."""
	observable_name = "Topological Charge"
	observable_name_compact = "topc"
	y_label = r"$Q$"
	x_label = r"$\sqrt{8t}$[fm]"
	formula = r"$Q = - \sum_x \frac{1}{64 \cdot 32\pi^2}\epsilon_{\mu\nu\rho\sigma}Tr\{G^{clov}_{\mu\nu}G^{clov}_{\rho\sigma}\}$"

class TopCharge4PostAnalysis(_PostAnalysis):
	"""Post-analysis of the topsus with Q^4."""
	observable_name = r"$\chi(\langle Q^4 \rangle)^{1/8}$"
	observable_name_compact = "topq4"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi(\langle Q^4 \rangle)^{1/8} [GeV]$" # 1/8 correct?
	formula = r"$\chi(\langle Q^4 \rangle)^{1/8} = \frac{\hbar}{aV^{1/4}} \langle Q^4 \rangle^{1/8} [GeV]$" # 1/8 correct?

class TopChargeTPostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topological charge at fixed euclidean time."""
	observable_name = "Topological Charge in Euclidean Time"
	observable_name_compact = "topct"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q_{t_{euclidean}}$"
	sub_obs = True

class TopChargeEuclSplitPostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topological charge in euclidean time intervals."""
	observable_name = "Topological Charge in Euclidean Time"
	observable_name_compact = "topcte"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q$"
	sub_obs = True

class TopChargeMCSplitPostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topological charge in MC time intervals."""
	observable_name = "Topological Charge in MC Time"
	observable_name_compact = "topcMC"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$Q$"
	sub_obs = True

#### Topological susceptibility definitions ####
class QtQ0PostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topsus at a fixed flow time."""
	observable_name = r"$\chi(\langle Q_t Q_{t_0} \rangle)^{1/4}$"
	observable_name_compact = "qtqzero"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi(\langle Q_{t} Q_{t_0} \rangle)^{1/4} [GeV]$" # $\chi_t^{1/4}[GeV]$
	sub_obs = True

class TopSusTPostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topsus with with one Q at fixed euclidean time."""
	observable_name = "Topological Susceptibility in Euclidean Time"
	observable_name_compact = "topsust"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi(\langle Q_t Q_{t_{euclidean}} \rangle)^{1/4} [GeV]$"
	sub_obs = True

class TopSusEuclSplitPostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topsus in euclidean time intervals."""
	observable_name = "Topological Susceptibility in MC Time"
	observable_name_compact = "topsuste"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi^{1/4} [GeV]$"
	sub_obs = True

class TopSusMCSplitPostAnalysis(_MultiObservablePostAnalysis):
	"""Post-analysis of the topsus in MC time intervals."""
	observable_name = "Topological Susceptibility in MC Time"
	observable_name_compact = "topsusMC"
	x_label = r"$\sqrt{8t_{flow}}[fm]$"
	y_label = r"$\chi^{1/4} [GeV]$"
	sub_obs = True


def main(args):
	"""
	Args should be post-analysis folder
	"""

	# Loads data from post analysis folder
	data = PostAnalysisDataReader(args[0])
	print "Retrieving data from folder: %s" % args[0]

	# Plots topsus
	topsus_analysis = TopSusPostAnalysis(data)
	topsus_analysis.set_analysis_data_type("bootstrap")
	topsus_analysis.plot()

	# Retrofits the topsus for continuum limit
	continium_targets = [0.3, 0.4, 0.5, 0.58]
	for cont_target in continium_targets:
		topsus_analysis.plot_continuum(cont_target)

	# Plots energy
	energy_analysis = EnergyPostAnalysis(data)
	energy_analysis.set_analysis_data_type("bootstrap")
	energy_analysis.plot()

	# Retrofits the energy for continuum limit
	energy_analysis.plot_continuum(0.3, 0.015, "bootstrap_fit")

	# Plot running coupling
	energy_analysis.coupling_fit()


	# PLAN
	# 1. Specify paths to data based data batch and batch name
	# 2. Load data in paths
	# 3. Specify what to plot. E.g. {observable} {analysis type} (optional) {euclidean or mc time}
	# 4. Plot




if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		args = sys.argv[1:]
	else:
		args = ["data5"]
		# args = ["dataGiovanni"]
		# args = ["../output/post_analysis_data/data4"]
	main(args)