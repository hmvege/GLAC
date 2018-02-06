"""
Functions that can be pickled by the multiprocessing module
"""
from statistics.autocorrelation import Autocorrelation
from statistics.jackknife import Jackknife
from statistics.bootstrap import Bootstrap


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
