from flow_analyser import *
from tools.folderreadingtools import DataReader
import statistics.parallel_tools as ptools
import os
import numpy as np
import copy
import sys
import time

def analyse_default(analysis_object, N_bs):
	analysis_object.boot(N_bs)
	analysis_object.jackknife()
	analysis_object.plot_original()
	analysis_object.plot_boot()
	analysis_object.plot_jackknife()
	analysis_object.autocorrelation()
	analysis_object.plot_autocorrelation(0)
	analysis_object.plot_autocorrelation(-1)
	analysis_object.plot_mc_history(0)
	analysis_object.plot_mc_history(-1)
	analysis_object.plot_original()
	analysis_object.plot_boot()
	analysis_object.plot_jackknife()
	analysis_object.plot_histogram(0)
	analysis_object.plot_histogram(-1)
	analysis_object.plot_integrated_correlation_time()
	analysis_object.plot_integrated_correlation_time()
	analysis_object.save_post_analysis_data()

def analyse_plaq(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params 
	plaq_analysis = AnalysePlaquette(obs_data("plaq"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
	analyse_default(plaq_analysis, N_bs)

def analyse_energy(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params 
	energy_analysis = AnalyseEnergy(obs_data("energy"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
	analyse_default(energy_analysis, N_bs)

def analyse_topsus(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params 
	topsus_analysis = AnalyseTopologicalSusceptibility(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
	analyse_default(topsus_analysis, N_bs)

def analyse_topc(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params 
	topc_analysis = AnalyseTopologicalCharge(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)

	if topc_analysis.beta == 6.0:
		topc_analysis.y_limits = [-9, 9]
	elif topc_analysis.beta == 6.1:
		topc_analysis.y_limits = [-12, 12]
	elif topc_analysis.beta == 6.2:
		topc_analysis.y_limits = [-12, 12]
	else:
		topc_analysis.y_limits = [None, None]

	analyse_default(topc_analysis, N_bs)

def analyse_topct(params):
	pass

def analyse_topc4(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params
	topcq4_analysis = AnalyseQQuartic(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
	analyse_default(topcq4_analysis, N_bs)

def analyse_qtq0(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params
	qzero_flow_time = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99] # Percents of data where we do qtq0
	qtqzero_analysis = AnalyseQtQZero(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)
	for qzero_flow_time_index in qzero_flow_time:
		# qtqzero_analysis.y_limits = [0, 2]
		qtqzero_analysis.setQ0(qzero_flow_time_index, y_label=r"$\langle Q_{t}Q_{t_0} \rangle^{1/4} [GeV]$")

		qtqzero_analysis.boot(N_bs)
		qtqzero_analysis.jackknife()
		qtqzero_analysis.plot_original()
		qtqzero_analysis.plot_jackknife()
		qtqzero_analysis.plot_boot()


def main(args):
	batch_name = args[0]
	batch_folder = args[1]
	obs_data = DataReader(batch_name, batch_folder, load_file=None, exclude_fobs=["topct"], dryrun=dryrun, verbose=verbose)

	obs_file = os.path.join("..",batch_folder,batch_name)
	os.path.isfile()

	N_bs = 500
	dryrun = False
	verbose = True
	parallel = True
	numprocs = 8
	paramers = [obs_data, dryrun, parallel, numprocs, verbose, N_bs]

	# fpath = obs_data.write_single_file()
	# obs_data = DataReader(directory_tree)
	# obs_data.retrieve_observable_data(exclude_fobs=["topct"])

	# # fpath = "../DataGiovanni/24_6.00.npy"
	# # fpath = "../DataGiovanni/28_6.10.npy"
	# fpath = "../DataGiovanni/32_6.20.npy"
	# NCfgs = 500
	# batch_folder = "DataGiovanni"
	# batch_name = "beta62"
	# # fpath = "../data5/beta60/24_6.00.npy"
	# obs_data = DataReader(batch_name, batch_folder, load_file=fpath, NCfgs=NCfgs, exclude_fobs=["topct"], dryrun=dryrun, verbose=verbose)

	# Analysis timers
	pre_time = time.clock()
	observable_strings = []

	post_time = time.clock()
	print "="*100
	print "Analysis of batch %s observables %s in %.2f seconds" % (batch_name, ", ".join([i.lower() for i in args[2:]]), (post_time-pre_time))
	print "="*100

if __name__ == '__main__':
	#### Available observables
	all_observables = ["plaq", "energy", "topc", "topsus", "qtqzero", "topc4", "topcq4", "topcqtq0"]
	basic_observables = ["plaq", "energy", "topc", "topsus"]

	#### Try to load binary file(much much faster)
	load_file = True

	#### Save binary file
	save_to_binary = True

	#### Different batches
	# data_batch_folder = "data2"
	# data_batch_folder = "data4"
	data_batch_folder = "data5"
	# giovanni_batch_folder = "DataGiovanni" 

	#### Different beta values folders:
	beta_folder_name = ["beta60", "beta61", "beta62"]
	# beta_folder_name = ["beta_60", "beta_61", "beta_62"]

	#### Basic batch setup
	databeta60 = {"batch_name": beta_folder_name[0],"batch_folder": data_batch_folder, 
		"observables": basic_observables, "NCfgs": 1000, "obs_file": "24_6.00",
		"load_file": load_file, "save_to_binary": save_to_binary}
	databeta61 = {"batch_name": beta_folder_name[1], "batch_folder": data_batch_folder, 
		"observables": basic_observables, "NCfgs": 500, "obs_file": "28_6.10",
		"load_file": load_file, "save_to_binary": save_to_binary}
	databeta62 = {"batch_name": beta_folder_name[2], "batch_folder": data_batch_folder,
		"observables": basic_observables, "NCfgs": 500, "obs_file": "32_6.20",
		"load_file": load_file, "save_to_binary": save_to_binary}

	#### Adding relevant batches to args
	args = [databeta60, databeta61, databeta62]

	#### Submitting observable-batches
	for a in args:
		main(a)
