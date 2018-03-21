import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import types

class LineFit:
	"""
	Line fit based on article by Lavagini et al (2007).
	"""
	def __init__(self, x, y, y_err=None):
		assert len(x) == len(y), "array lengths to not match"
		self.x = x
		self.y = y
		self.y_err = y_err
		self.n = len(y)
		self.x_lower = np.min(self.x)
		self.x_upper = np.max(self.x)

	def fit(self, x_arr=None):
		"""
		Fits a linear function,
			y_hat = b0 + b1 * x
		from a provided dataset, x, y and y_err from initialization.

		Args:
			x_arr: optional array for which return as a fit. Default is None.

		Returns:
			y_hat: returns fit(x)
			y_hat_err: returns fit_error(x)
			fit_params: list containing b0, b0_err, b1, b1_err
			chi_squared: goodness of fit
		"""

		if isinstance(x_arr, types.NoneType):
			x_arr = np.linspace(self.x[self.x_lower], self.x[self.x_upper], 100)

		self.x_mean = np.mean(self.x)
		self.y_mean = np.mean(self.y)

		# Temporary som contained in both eq. 4, eq. 6 and eq. 7
		self.xi_xmean_sum = np.sum((self.x - self.x_mean)**2)

		# Eq. 4
		self.b1 = np.sum((self.x - self.x_mean) * self.y) / self.xi_xmean_sum

		# Eq. 3
		self.b0 = self.y_mean - self.b1 * self.x_mean

		# Eq. 5
		self.s_xy_err = np.sum((self.y - self._y_hat(self.x))**2) / (self.n - 2)

		# Eq. 6
		b0_err = self.s_xy_err
		b0_err *= (1.0 / self.n + self.x_mean ** 2 / self.xi_xmean_sum)

		# Eq. 7
		b1_err = self.s_xy_err / self.xi_xmean_sum

		fit_params = [self.b0, b0_err, self.b1, b1_err]

		# Goodness of fit
		chi = self._chi_squared(self.y, self.y_err, self._y_hat(self.x))

		return self._y_hat(x_arr), self._y_hat_err(x_arr), fit_params, chi

	def _y_hat(self, x):		
		"""Unweigthed y(x), eq. 1."""
		return self.b0 + self.b1 * x

	def _y_hat_err(self, x):
		"""Unweigthed y(x) error, eq. 8b."""
		_pt1 = self.b0 + self.b1 * x
		_pt2 = scipy.stats.t.isf(0.68, self.n - 2) * np.sqrt(self.s_xy_err) 
		_pt2 *= np.sqrt(1.0 / self.n + (x - self.x_mean)**2 / self.xi_xmean_sum)
		return [_pt1 + _pt2, _pt1 - _pt2]

	def fit_weighted(self, x_arr=None):
		"""
		Performs a weighted fit based on x, y, and y_err provided in 
		initialization.

		Args:
			x_arr: optional array for which return as a fit. Default is None.

		Returns:
			y_hat: returns fit(x)
			y_hat_err: returns fit_error(x)
			fit_params: list containing b0, b0_err, b1, b1_err
			chi_squared: goodness of fit
		"""

		if isinstance(x_arr, types.NoneType):
			x_arr = np.linspace(self.x[0], self.x[-1], 100)

		assert not isinstance(self.y_err, types.NoneType), "Missing y_err."

		# w = self.n * 1.0/self.y_err**2 * 1.0/np.sum(self.y_err**2)
		self.w = 1.0/self.y_err**2

		self.xw_mean = np.sum(self.w * self.x) / np.sum(self.w)
		self.yw_mean = np.sum(self.w * self.y) / np.sum(self.w)

		# Temporary som contained in both eq. 4, eq. 6 and eq. 7
		self.xwi_xmean_sum = np.sum(self.w * (self.x - self.xw_mean)**2)

		# Eq. 18
		self.b1w = np.sum(self.w * (self.x - self.xw_mean) * self.y)
		self.b1w /= self.xwi_xmean_sum

		# Eq. 17
		self.b0w = self.yw_mean - self.b1w * self.xw_mean

		# Eq. 21
		self.s_xyw_err = np.sum(self.w * (self.y - self._yw_hat(self.x))**2)
		self.s_xyw_err /= (self.n - 2.0) 

		# Eq. 19
		b0w_err = (1.0/np.sum(self.w) + self.xw_mean**2 / self.xwi_xmean_sum) 
		b0w_err *= self.s_xyw_err

		# Eq. 20
		b1w_err = self.s_xyw_err / self.xwi_xmean_sum

		fit_params = [self.b0w, b0w_err, self.b1w, b1w_err]

		# Goodness of fit
		chi = self._chi_squared(self.y, self.y_err, self._yw_hat(self.x))

		return self._yw_hat(x_arr), self._yw_hat_err(x_arr), fit_params, chi

	def _yw_hat(self, x):
		"""weigthed y(x), eq. 1"""
		return self.b0w + self.b1w * x

	def _yw_hat_err(self, x):
		"""Weigthed y(x) errors, eq. 22."""
		_pt1 = self.b0w + self.b1w * x
		_pt2 = scipy.stats.t.isf(0.68, self.n - 2) * np.sqrt(self.s_xyw_err) 
		_pt2 *= np.sqrt(1.0 / np.sum(self.w) + (x - self.xw_mean)**2 / self.xwi_xmean_sum)
		return [_pt1 + _pt2, _pt1 - _pt2]

	def inverse_fit(self, y0, weigthed=False):
		"""
		Inverse fiting.

		Args:
			y0: target fit at y-axis, float.
			weigthed: bool, if we are to use weighted fit or not.

		Returns:
			x0: targeted y0 fit
			x0_err: errorband at y0 fit
		"""
		_n = 10000
		_x = np.linspace(self.x_lower, self.x_upper, _n)

		if weigthed:
			# Finds the target value
			x0_index = np.argmin(np.abs(y0 - self._yw_hat(_x)))
			x0 = _x[x0_index]

			# Finds the error bands
			_x_err_pos, _x_err_neg = self._yw_hat_err(_x)
			x0_err_index_neg = np.argmin(np.abs(y0 - _x_err_pos))
			x0_err_index_pos = np.argmin(np.abs(y0 - _x_err_neg))
			x0_err = [_x[x0_err_index_neg], _x[x0_err_index_pos]]

		else:
			# Finds the target value
			min_index = np.argmin(np.abs(y0 - self._y_hat(_x)))
			x0 = _x[min_index]

			# Finds the error bands
			_x_err_pos, _x_err_neg = self._y_hat_err(_x)
			x0_err_index_neg = np.argmin(np.abs(y0 - _x_err_pos))
			x0_err_index_pos = np.argmin(np.abs(y0 - _x_err_neg))
			x0_err = [_x[x0_err_index_neg], _x[x0_err_index_pos]]

		return x0, x0_err

	@staticmethod
	def _chi_squared(y, y_err, y_fit):
		"""Goodness test of fit."""
		return np.sum((y - y_fit)**2 / y_err**2)

def main():
	import random
	import scipy.optimize as sciopt

	def fit_var_printer(var_name, var, var_error, w=16):
		return "{0:<s} = {1:<.{w}f} +/- {2:<.{w}f}".format(var_name, var, var_error, w=w)

	# Generates signal with noise
	a = 0.65
	b = 1.1
	N = 4
	x = np.linspace(0, 5, N)
	signal_spread = 0.5
	signal = a*x + b + np.random.uniform(-signal_spread, signal_spread, N)
	signal_std = np.random.uniform(0.1, 0.3, N)

	# x = np.append(x, 2.5)
	# signal = np.append(signal_std, 3.5)
	# signal_std = np.append(signal_std, 5.0)

	fit = LineFit(x, signal, signal_std)
	x_hat = np.linspace(0, 5, 100)

	# Unweigthed fit
	y_hat, y_hat_err, f_params, chi_unweigthed = fit.fit(x_hat)
	b0, b0_err, b1, b1_err = f_params

	# Fits without any weights first
	pol1, polcov1 = sciopt.curve_fit(lambda x, a, b : x*a + b, x, signal)

	# Fit target
	fit_target = 2.5
	x_fit, x_fit_error = fit.inverse_fit(fit_target)

	print "UNWEIGTHED LINE FIT"
	print "SciPy curve_fit:"
	print fit_var_printer("a", pol1[0], polcov1[0,0])
	print fit_var_printer("b", pol1[1], polcov1[1,1])


	print "LineFit:"
	print fit_var_printer("a", b1, b1_err)
	print fit_var_printer("b", b0, b0_err)
	print "Goodness of fit: %f" % chi_unweigthed
	# print "b = {0:<.10f} +/- {1:<.10f}".format(b0, b0_err)

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(211)
	ax1.axhline(fit_target, linestyle="dashed", color="tab:grey")
	ax1.plot(x_hat, y_hat, label="Unweigthed fit", color="tab:blue")
	ax1.fill_between(x_hat, y_hat_err[0], y_hat_err[1], alpha=0.5, color="tab:blue")
	ax1.plot(x, signal, marker="o", label="Signal", linestyle="none", color="tab:orange")

	ax1.set_ylim(0.9, 4.5)
	ax1.axvline(x_fit, color="tab:orange")
	ax1.fill_betweenx(np.linspace(0,6,100), x_fit_error[0], x_fit_error[1], label=r"$x_0\pm\sigma_{x_0}$", alpha=0.5, color="tab:orange")
	ax1.legend(loc="best", prop={"size":8})
	ax1.set_title("Fit test - unweigthed")

	# Weighted fit
	yw_hat, yw_hat_err, f_params_weigthed, chi_weigthed = fit.fit_weighted(x_hat)
	b0, b0_err, b1, b1_err = f_params_weigthed
	xw_fit, xw_fit_error = fit.inverse_fit(fit_target, weigthed=True)

	print "WEIGTHED LINE FIT"
	print fit_var_printer("a", b1, b1_err)
	print fit_var_printer("b", b0, b0_err)
	print "Goodness of fit: %f" % chi_weigthed

	ax2 = fig1.add_subplot(212)
	ax2.axhline(fit_target, linestyle="dashed", color="tab:grey")
	ax2.errorbar(x, signal, yerr=signal_std, fmt="o", label="Signal", color="tab:orange")
	ax2.plot(x_hat, yw_hat, label="Weighted fit", color="tab:blue")
	ax2.fill_between(x_hat, yw_hat_err[0], yw_hat_err[1], alpha=0.5, color="tab:blue")

	ax2.set_ylim(0.6, 5)
	ax2.axvline(xw_fit, color="tab:orange")
	ax2.fill_betweenx(np.linspace(0,6,100), xw_fit_error[0], xw_fit_error[1], label=r"$x_{0,w}\pm\sigma_{x_0,w}$", alpha=0.5, color="tab:orange")

	ax2.legend(loc="best", prop={"size":8})
	plt.show()

if __name__ == '__main__':
	main()