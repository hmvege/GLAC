import numpy as np, matplotlib.pyplot as plt, sys, os, re, scipy.optimize as sciopt

# from scipy.optimize import curve_fit

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a # fermi

# def linear_curve(x,a,b):
# 	return a*x + b

class LineFit:
	"""
	Line fit method when y contains errors
	"""
	def __init__(self, x, y, y_err):
		None

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

		for value in plot_values:
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
		# plt.show()
		plt.close(fig)

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

		self.plot_values = sorted(plot_values,key=lambda x:x["batch_name"])

		# x_limits = [0,np.max(plot_values[-1]["x"])]
		# y_limits = [0,np.max(plot_values[-1]["y"])]

		self._plot_core(self.plot_values)#,x_limits=x_limits,y_limits=y_limits)

	def _get_t0(self,data_values):
		self.t0_values = []
		for values in data_values:
			t0_batch = {}
			t0_batch["beta"] = self.get_float(values["batch_name"])
			t0_batch["t0"],t0_batch["t0_err"] = self._linefit_t0(values["x"],values["y"],values["y_err"],t0_batch["beta"])
			t0_batch["a"] = getLatticeSpacing(t0_batch["beta"])

			# Adds to list of batch-values
			self.t0_values.append(t0_batch)

	def _linefit_t0(self,x,y,y_err,beta):
		# Creates a small interval in which we will fit the data
		y0 = 0.3
		fit_interval = 0.015

		x_line_to_fit = []
		y_line_to_fit = []
		y_line_to_fit_errors = []
		for i, iy in enumerate(y):
			if (y0 - fit_interval) < iy < (y0 + fit_interval):
				x_line_to_fit.append(x[i])
				y_line_to_fit.append(iy)
				y_line_to_fit_errors.append(y_err[i])

		if len(y_line_to_fit) == 0:
			raise ValueError("No values in the vincinity of t^2<E>=0.3 found for beta %2.2f" % beta)

		x_line_to_fit = np.asarray(x_line_to_fit)
		y_line_to_fit = np.asarray(y_line_to_fit)
		y_line_to_fit_errors = np.asarray(y_line_to_fit_errors)
		# print x_line_to_fit,y_line_to_fit,y_line_to_fit_errors

		# Fitting data
		pol,polcov = np.polyfit(x_line_to_fit,y_line_to_fit,1,rcond=None,full=False,w=1.0/y_line_to_fit_errors,cov=True)
		# print pol,"\n",polcov
		# pol, polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x_line_to_fit, y_line_to_fit, sigma = y_line_to_fit_errors)
		# print pol,"\n",polcov

		# Retrieving polynomial values for retrofitting
		a = pol[0]
		b = pol[1]
		a_err, b_err = np.sqrt(np.diag(polcov))

		# t0 value
		t0 = self.r0**2 * (y0 - b) / a

		# Gets error of t0 estimate
		t0_err = np.sqrt((self.r0**2*(b - y0)/a**2)**2 * a_err**2 + self.r0**4/a**2 * b_err**2 + 2*self.r0**4*(b - y0)/a**4 * polcov[1,0])

		# print "beta %g:" % beta, t0_err/self.r0**2
		# print "beta %g:" % beta, t0/self.r0**2, y0

		# plt.errorbar(x_line_to_fit,y_line_to_fit,yerr=y_line_to_fit_errors)
		# plt.errorbar(t0/self.r0**2,y0,yerr=t0_err/self.r0**2,color="r",ecolor="r")
		# plt.show()
		# exit(1)

		return t0,t0_err

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
		# pol, polcov = np.polyfit(x_datapoints,y_datapoints,1,rcond=None,full=False,w=1.0/y_datapoints_error,cov=True)
		# print pol,"\n",polcov

		pol, polcov = sciopt.curve_fit(lambda x, a, b : x*a + b, x_datapoints, y_datapoints,sigma=y_datapoints_error)
		# print pol,"\n",polcov

		# exit(1)

		# Gets line properties
		a = pol[0]
		b = pol[1]
		a_err, b_err = np.sqrt(np.diag(polcov))

		# Sets up line fitted variables		
		x = np.linspace(0,x_datapoints[-1]*1.03,100)
		y = a * x + b
		y_std = a_err * x + b_err

		return x,y,y_std,x_datapoints,y_datapoints,y_datapoints_error, a, b, a_err, b_err

	def plot_continiuum(self):
		# Retrieves t0 values
		self._get_t0(self.plot_values)

		# Fits to continiuum and retrieves values to be plotted
		x, y, y_std, ar0, t0, t0_err, a, b, a_err, b_err = self._linefit_to_continiuum()

		x_points = np.zeros(len(ar0)+1)
		y_points = np.zeros(len(ar0)+1)
		y_points_err = np.zeros(len(ar0)+1)

		x_points[0] = x[0]
		y_points[0] = y[0]
		y_points_err[0] = y_std[0]

		x_points[1:] = ar0
		y_points[1:] = t0
		y_points_err[1:] = t0_err

		fig = plt.figure(self.dpi)
		ax = fig.add_subplot(111)

		# print x_points,y_points,y_points_err

		# ax.axvline(0,linestyle="--",color="0",alpha=0.5)
		ax.errorbar(x_points[1:],y_points[1:],yerr=y_points_err[1:],fmt="o",color="0",ecolor="0")
		ax.errorbar(x_points[0],y_points[0],yerr=y_points_err[0],fmt="o",capthick=4,color="r",ecolor="r")
		ax.plot(x,y,color="0")#,label=r"$y=%2.4fx + %2.4f$" % (a,b))
		ax.set_ylabel(r"$\frac{\sqrt{8t_0}}{r_0}$")
		ax.set_xlabel(r"$(a/r_0)^2$")
		# ax.set_title(r"Continiuum limit reference scale: $t_{0,cont}=%2.4f\pm%g$" % ((self.r0*y_points[0])**2/8,(self.r0*y_points_err[0])**2/8))
		ax.set_xlim(-0.005,0.045)
		ax.set_ylim(0.92,0.98)


		start, end = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(start, end, 0.01))
		ax.grid(True)

		ax.annotate(r"$a=0.05$fm", xy=((0.01/self.r0)**2, end), xytext=((0.01/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")
		ax.annotate(r"$a=0.07$fm", xy=((0.07/self.r0)**2, end), xytext=((0.07/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")
		ax.annotate(r"$a=0.07$fm", xy=((0.1/self.r0)**2, end), xytext=((0.1/self.r0)**2, end+0.005),arrowprops=dict(arrowstyle="->"),ha="center")

		# ax.legend(loc="lower left")

		fname = os.path.join(self.output_folder,"post_analysis_%s_continiuum.png" % self.analysis_name_compact)
		fig.savefig(fname,dpi=self.dpi)

		print "Figure created in %s" % fname
		plt.show()
		plt.close()


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
	energy_analysis.plot_continiuum()
	

if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		args = sys.argv[1:]
	else:
		args = ["../output/post_analysis_data"]
	main(args)