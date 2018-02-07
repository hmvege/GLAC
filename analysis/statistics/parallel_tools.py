from autocorrelation import Autocorrelation
from jackknife import Jackknife
from bootstrap import Bootstrap
import numpy as np

"""
Functions that can be pickled by the multiprocessing module
"""

def _autocorrelation_parallel_core(input_values):
	ac = Autocorrelation(input_values[0], use_numpy = input_values[1], data_statistic = input_values[2])
	return ac.R, ac.R_error, ac.integrated_autocorrelation_time(), ac.integrated_autocorrelation_time_error()

def _bootstrap_parallel_core(input_values):
	data, N_bs, bs_statistic, non_bs_stats, index_lists = input_values
	bs = Bootstrap(data, N_bs, bootstrap_statistics = bs_statistic, non_bs_stats = non_bs_stats, index_lists = index_lists)
	return bs.bs_avg, bs.bs_std, bs.avg_original, bs.std_original, bs.bs_data, bs.data_original

def _jackknife_parallel_core(input_values):
	data, jk_statistics, non_jk_statistics = input_values
	jk = Jackknife(data, jk_statistics = jk_statistics, non_jk_statistics = non_jk_statistics)
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

def _chi_beta6_0(Q_squared):
	const = 0.001152302372957481
	return const*Q_squared**(0.25)

def _chi_beta6_0_error(Q_squared, Q_squared_std):
	const = 0.001152302372957481
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_1(q_squared):
	const = 0.0009714961458189535
	return const*Q_squared**(0.25)

def _chi_beta6_1_error(Q_squared, Q_squared_std):
	const = 0.001152302372957481
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_2(q_squared):
	const = 0.0008363484965013975
	return const*Q_squared**(0.25)

def _chi_beta6_2_error(Q_squared, Q_squared_std):
	const = 0.001152302372957481
	return 0.25*const*Q_squared_std / Q_squared**(0.75)

def _chi_beta6_45(q_squared):
	const = 0.0005359545920732469
	return const*Q_squared**(0.25)

def _chi_beta6_45_error(Q_squared, Q_squared_std):
	const = 0.001152302372957481
	return 0.25*const*Q_squared_std / Q_squared**(0.75)