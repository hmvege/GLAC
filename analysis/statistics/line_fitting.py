import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt

def fit_line_form_bootstrap(x, bs_data, observable, beta, fit_target, fit_interval, axis = "y", fit_function_modifier = lambda x : x, fit_method = "curve_fit", plot_fit_window = False):
	"""
	Function for creating a line fit using bootstrapped data.
	Args:
		x 					(numpy float array): x-axis values to fit along
		bs_data					   (dictionary): array of un-averaged bootstrapped data. NFlows x NBootstraps
		observable 						  (str): name of the observable
		beta 							(float): beta value
		axis 							  (str): axis to fit for. Options = "x","y"(default)
		fit_target						(float): y-axis value to fit against
		fit_interval					(float): y-axis interval in which we will linefit for
		[fit_function_modifier] 	 (function): function for modifying y-values
		[fit_method]					  (str): type of data fit to be performed. Default is scipy.optimize.curve_fit
		[plot_fit_window]				 (bool): plots the fit window to easier visualize
	Returns:
		fit_value 						(float): y value fitted at fit_target
		fit_value_error 				(float): error of y_0
	"""
	# Retrieves data
	data = bs_data
	NFlows, NBoot = data.shape
	# print data.shape

	# # Applies the funciton correction to all bootstrap samples
	# for iBoot in xrange(NBoot):
	# 	data[:,iBoot] = data_function_correction(data[:,iBoot])

	# Sets the fitting method for the each bootstrap value
	fit_methods_list = ["curve_fit","polynomial"]
	if fit_method not in fit_methods_list:
		raise KeyError("%s not a possible fit method. Use %s" % (fit_method,", ".join(fit_methods_list)))

	# Takes data mean for finding fit interval
	data_mean = np.mean(bs_data,axis=1)
	data_std = np.std(bs_data,axis=1)

	if axis == "y":
		# Finds fit interval indices
		start_index = np.argmin( np.abs(data_mean - (fit_target - fit_interval)) )
		end_index = np.argmin( np.abs(data_mean - (fit_target + fit_interval)) )
	elif axis == "x":
		# Finds fit interval indices against y axis values
		start_index = np.argmin( np.abs(x - (fit_target - fit_interval)) )
		end_index = np.argmin( np.abs(x - (fit_target + fit_interval)) )
	else:
		raise KeyError("%s is not a valid axis. Options: 'x','y'(default)" % axis)

	N_fit_line_length = end_index - start_index

	# Plots the fit window if prompted to do so
	if plot_fit_window:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(x,data_mean,label=observable)
		ax.fill_between(x, data_mean - data_std, data_mean + data_std,alpha=0.5)
		ax.scatter(x[start_index],data_mean[start_index],label="Fit line start")
		ax.scatter(x[end_index],data_mean[end_index],label="Fit line end")
		ax.axhline(fit_target,color="0",alpha=0.5)
		ax.grid(True)
		ax.legend()
		ax.set_title(r"Fit window for %s at $\beta=%.2f$" % (observable,beta))
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

		if axis == "y":
			# Adds bootstrap fit to fitted values if we are fitting y axis
			fitted_values[iBoot] = fit_function_modifier((fit_target - b) / a)
		else:
			# Adds bootstrap fit to fitted values if we are fitting for x axis
			fitted_values[iBoot] = fit_function_modifier(fit_target*a + b)

	# Takes mean and standard deviation of fitted values
	fit_value = np.mean(fitted_values)
	fit_value_error = np.std(fitted_values)

	return fit_value, fit_value_error

def fit_line(x, y, y_error, observable, beta, fit_target, fit_interval, axis="y", fit_function_modifier=lambda x: x, fit_method="curve_fit", plot_fit_window=False):
	"""
	Function for creating a line fit to a fit target and within a fit interval.
	Args:
		x 					(numpy float array): x values to fit along
		y 					(numpy float array): y values to fit along. length of NFlows
		y_error				(numpy float array): y value errors
		observable 						  (str): name of the observable
		beta 							(float): beta value
		axis 							  (str): axis to fit for. Options = "x","y"(default)
		fit_target						(float): y-axis value to fit against
		fit_interval					(float): y-axis interval in which we will linefit for
		[fit_function_modifier] 	 (function): function for modifying y-values
		[fit_method]					  (str): type of data fit to be performed. Default is scipy.optimize.curve_fit
		[plot_fit_window]				 (bool): plots the fit window to easier visualize
	Returns:
		fit_value 						(float): y value fitted at fit_target
		fit_value_error 				(float): error of y_0
	"""
	NFlows = len(y)

	# Small error catching ensuring equal length of arrays
	assert len(x) == len(y) == len(y_error), "arrays x, y and y_error do not match."

	# Sets the fitting method for the each bootstrap value
	fit_methods_list = ["curve_fit","polynomial"]
	if fit_method not in fit_methods_list:
		raise KeyError("%s not a possible fit method. Use %s" % (fit_method,", ".join(fit_methods_list)))

	if axis == "y":
		# Finds fit interval indices against y axis values
		start_index = np.argmin( np.abs(y - (fit_target - fit_interval)) )
		end_index = np.argmin( np.abs(y - (fit_target + fit_interval)) )
	elif axis == "x":
		# Finds fit interval indices against y axis values
		start_index = np.argmin( np.abs(x - (fit_target - fit_interval)) )
		end_index = np.argmin( np.abs(x - (fit_target + fit_interval)) )
	else:
		raise KeyError("%s is not a valid axis. Options: 'x','y'(default)" % axis)

	N_fit_line_length = end_index - start_index

	# Plots the fit window if prompted to do so
	if plot_fit_window:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(x, y, label=observable)
		ax.fill_between(x, y - y_error, y + y_error, alpha=0.5)
		ax.scatter(x[start_index], y[start_index], label="Fit line start")
		ax.scatter(x[end_index], y[end_index], label="Fit line end")
		if axis == "y":
			ax.axhline(fit_target, color="0", alpha=0.5)
		else:
			ax.axvline(fit_target, color="0", alpha=0.5)
		ax.grid(True)
		ax.legend()
		ax.set_title(r"Fit window for %s at $\beta=%.2f$" % (observable, beta))
		plt.show()
		plt.close(fig)

	# Fitting data
	if fit_method == "curve_fit":
		pol,polcov = sciopt.curve_fit(lambda x, a, b: x*a + b,
			x[start_index:end_index], y[start_index:end_index],
			sigma=y_error[start_index:end_index])

	elif fit_method == "polynomial":
		pol, polcov = np.polyfit(x[start_index:end_index],
			y[start_index:end_index,iBoot], 1, rcond=None, full=False, 
			cov=True, w=1.0/(y_error[start_index:end_index]))
	else:
		raise KeyError("No fit method called %s found." % fit_method)

	# Retrieving polynomial values for retrofitting
	a = pol[0]
	b = pol[1]
	a_err, b_err = np.sqrt(np.diag(polcov))

	if axis == "y":
		# Gets the fit values
		fit_value = fit_function_modifier((fit_target - b)/a)
		fit_value_error = fit_function_modifier(fit_value*np.sqrt((b_err/(fit_target - b))**2 + (a_err/a)**2))
	else:
		# Adds bootstrap fit to fitted values if we are fitting for x axis
		fit_value = fit_function_modifier(fit_target*a + b)
		fit_value_error = fit_function_modifier(fit_value*a_err + b_err)

	return fit_value, fit_value_error

if __name__ == '__main__':
	exit("Exit: %s to be imported as module in other programs." % __file__.split("/")[-1])