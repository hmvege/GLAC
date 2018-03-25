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
	data, funder = input_values
	if funder == None:
		funder = lambda x: x
	ac = PropagatedAutocorrelation(data, function_derivative=funder)
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

def _return_mean_squared(x,axis=None):
	# For use instead of lambda x**2 : np.mean(x**2) in parallel
	return np.mean(x**2,axis=axis)

# const = hbarc / a(beta) / (V(beta)**(0.25))

# BETA 6.0
def _chi_beta6_0(Q_squared):
	const = 0.074230020266
	return const*Q_squared**(0.25)

def _chi_beta6_0_error(Q_squared, Q_squared_std):
	const = 0.074230020266
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_0_derivative(Q_squared):
	const = 0.074230020266
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_0_derivative_int12(Q_squared):
	# When we are looking at euclidean intervals of size 16
	const = 0.074230020266 / (0.25)**(0.25)
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_0_derivativeQ4(Q4):
	const = 0.074230020266
	return (0.25*const)**2 / Q4**(0.5625)

# BETA 6.1
def _chi_beta6_1(Q_squared):
	const = 0.0749573791735
	return const*Q_squared**(0.25)

def _chi_beta6_1_error(Q_squared, Q_squared_std):
	const = 0.0749573791735
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_1_derivative(Q_squared):
	const = 0.0749573791735
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_1_derivative_int14(Q_squared):
	# When we are looking at euclidean intervals of size 14
	const = 0.0749573791735 / (0.25)**(0.25)
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_1_derivativeQ4(Q4):
	const = 0.0749573791735
	return (0.25*const)**2 / Q4**(0.5625)

# BETA 6.2
def _chi_beta6_2(Q_squared):
	const = 0.0763234462734
	return const*Q_squared**(0.25)

def _chi_beta6_2_error(Q_squared, Q_squared_std):
	const = 0.0763234462734
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_2_derivative(Q_squared):
	const = 0.0763234462734
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_2_derivative_int16(Q_squared):
	# When we are looking at euclidean intervals of size 16
	const = 0.0763234462734 / (0.25)**(0.25)
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_2_derivativeQ4(Q4):
	const = 0.0763234462734
	return (0.25*const)**2 / Q4**(0.5625)

# BETA 6.45
def _chi_beta6_45(Q_squared):
	const = 0.0723048176484
	return const*Q_squared**(0.25)

def _chi_beta6_45_error(Q_squared, Q_squared_std):
	const = 0.0723048176484
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_45_derivative(Q_squared):
	const = 0.0723048176484
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_45_derivative_int24(Q_squared):
	# When we are looking at euclidean intervals of size 24
	const = 0.0723048176484 / (0.25)**(0.25)
	return 0.25*const / Q_squared**(0.75)

def _chi_beta6_45_derivativeQ4(Q4):
	const = 0.0723048176484
	return (0.25*const)**2 / Q4**(0.5625)


if __name__ == '__main__':
	exit("Exit: %s to be imported as module in other programs." % __file__.split("/")[-1])