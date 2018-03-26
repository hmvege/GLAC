from autocorrelation import Autocorrelation, PropagatedAutocorrelation
from jackknife import Jackknife
from bootstrap import Bootstrap
import numpy as np

"""
Functions that can be pickled by the multiprocessing module
"""

def _autocorrelation_parallel_core(input_values):
	ac = Autocorrelation(input_values[0])
	return ac.R, ac.R_error, ac.integrated_autocorrelation_time(), ac.integrated_autocorrelation_time_error()

def _autocorrelation_propagated_parallel_core(input_values):
	data, funder, funder_params = input_values
	if funder == None:
		funder = lambda x: x
	assert isinstance(funder_params, dict), "function parameters is not a dictionary."
	ac = PropagatedAutocorrelation(data, function_derivative=funder, function_parameters=funder_params)
	return ac.R, ac.R_error, ac.integrated_autocorrelation_time(), ac.integrated_autocorrelation_time_error()

def _bootstrap_parallel_core(input_values):
	data, N_bs, index_lists = input_values
	bs = Bootstrap(data, N_bs, index_lists=index_lists)
	return bs.bs_avg, bs.bs_std, bs.avg_original, bs.std_original, bs.bs_data, bs.data_original

def _jackknife_parallel_core(input_values):
	jk = Jackknife(input_values)
	return jk.jk_avg, jk.jk_std, jk.jk_data

def _default_return(x):
	# For use instead of lambda x : x in parallel
	return x

def _default_error_return(x, x_std):
	# For use instead of lambda x : x in parallel
	return x_std

def _return_squared(x):
	# For use instead of lambda x**2 : x**2 in parallel
	return x*x

def _return_mean_squared(x, axis=None):
	# For use instead of lambda x**2 : np.mean(x**2) in parallel
	return np.mean(x**2, axis=axis)

# Functions to be passed with **kwargs
def _chi(Q2, const=None):
	return const*Q2**(0.25)

def _chi_error(Q2, Q2_std, const=None):
	return 0.25*const*Q2_std / Q2**(0.75)

def _chi_derivative(Q2, const=None):
	return 0.25*const / Q2**(0.75)

def _chi_derivativeQ4(Q4, const=None):
	return (0.25*const)**2 / Q4**(0.5625)

if __name__ == '__main__':
	exit("Exit: %s to be imported as module in other programs." % __file__.split("/")[-1])