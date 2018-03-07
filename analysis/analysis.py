from flow_analyser import *
from tools.folderreadingtools import DataReader
import statistics.parallel_tools as ptools
import os
import numpy as np
import copy
import sys
import time

def analyse_default(analysis_object, N_bs):
	print analysis_object
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

def analyse_topc4(params):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params
	topcq4_analysis = AnalyseQQuartic(obs_data("topc"), dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)

	analyse_default(topcq4_analysis, N_bs)

def analyse_qtq0(params, qzero_flow_times):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	qtqzero_analysis = AnalyseQtQZero(obs_data("topc"), dryrun=dryrun,
		parallel=parallel, numprocs=numprocs, verbose=verbose)

	for qzero_flow_time_index in qzero_flow_times:
		qtqzero_analysis.setQ0(qzero_flow_time_index)
		analyse_default(qtqzero_analysis, N_bs)

def analyse_topct(params, t_euclidean_indexes):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	topct_analysis = AnalyseTopologicalChargeInEuclideanTime(obs_data("topct"),
		dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)

	for ie in t_euclidean_indexes[topct_analysis.beta]:
		topct_analysis.setEQ0(ie)
		analyse_default(topct_analysis, N_bs)

def analyse(parameters):
	"""
	Function for starting flow analyses.

	Args:
		parameters: dictionary containing following elements: batch_name, 
			batch_folder, observables, NCfgs, obs_file, load_file, 
			save_to_binary, base_parameters, flow_epsilon, NFlows,
			create_perflow_data, correct_energy
	"""

	# Analysis timers
	pre_time = time.clock()
	observable_strings = []

	# Retrieves analysis parameters
	batch_name = parameters["batch_name"]
	batch_folder = parameters["batch_folder"]

	_bp = parameters["base_parameters"]

	# Retrieves data
	obs_data = DataReader(batch_name, batch_folder,
		load_file=parameters["load_file"],
		flow_epsilon=parameters["flow_epsilon"], NCfgs=parameters["NCfgs"],
		create_perflow_data=parameters["create_perflow_data"],
		dryrun=_bp["dryrun"], verbose=_bp["verbose"],
		correct_energy=parameters["correct_energy"])

	# Writes raw observable data to a single binary file
	if parameters["save_to_binary"] and not parameters["load_file"]:
		obs_data.write_single_file()

	# Builds parameters list to be passed to analyser
	params = [obs_data, _bp["dryrun"], _bp["parallel"], _bp["numprocs"], 
		_bp["verbose"], _bp["N_bs"]]

	# Runs through the different observables and analyses each one
	if "plaq" in parameters["observables"]:
		analyse_plaq(params)
	if "energy" in parameters["observables"]:
		analyse_energy(params)
	if "topc" in parameters["observables"]:
		analyse_topc(params)
	if "topsus" in parameters["observables"]:
		analyse_topsus(params)
	if "qtqzero" in parameters["observables"]:
		analyse_qtq0(params, parameters["qzero_flow_times"])
	if "topc4" in parameters["observables"]:
		analyse_topc4(params)
	if "topct" in parameters["observables"]:
		analyse_topct(params, parameters["topct_euclidean_indexes"])
	
	post_time = time.clock()
	print "="*100
	print "Analysis of batch %s observables %s in %.2f seconds" % (batch_name,
		", ".join([i.lower() for i in parameters["observables"]]), (post_time-pre_time))
	print "="*100

def main():
	#### Available observables
	all_observables = ["plaq", "energy", "topc", "topsus", "qtqzero", "topc4", "topct"]
	basic_observables = ["plaq", "energy", "topc", "topsus"]
	# observables = all_observables[6:7]
	observables = all_observables[6:7] + ["qtqzero", "topc4"]
	# observables = basic_observables
	print "Running for observables: %s" % ", ".join(observables)

	#### Base parameters
	N_bs = 500
	dryrun = False
	verbose = True
	parallel = True
	numprocs = 8
	base_parameters = {"N_bs": N_bs, "dryrun": dryrun, "verbose": verbose, 
		"parallel": parallel, "numprocs": numprocs}

	#### Try to load binary file(much much faster)
	load_file = True

	# If we are to create per-flow datasets as opposite to per-cfg datasets
	create_perflow_data = False

	#### Save binary file
	save_to_binary = True

	#### Load specific parameters
	NFlows = 1000
	flow_epsilon = 0.01

	#### Different batches
	# data_batch_folder = "data2"
	# data_batch_folder = "data4"
	data_batch_folder = "data5"
	# data_batch_folder = "DataGiovanni" 

	#### If we need to multiply
	if data_batch_folder == "DataGiovanni":
		correct_energy = False
		load_file = True
		save_to_binary = False
	else:
		correct_energy = True

	#### Different beta values folders:
	beta_folder_name = ["beta60", "beta61", "beta62"]
	# beta_folder_name = ["beta_60", "beta_61", "beta_62"]

	# Indexes to look at for topct
	t_euclidean_indexes = {
		6.0: [0, 11, 23, 35, 47],
		6.1: [0, 13, 27, 41, 55],
		6.2: [0, 15, 31, 47, 63],
		6.45: [0, 23, 47, 71, 95]}

	# Percents of data where we do qtq0
	qzero_flow_times = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]

	#### Basic batch setup
	default_params = {"batch_folder": data_batch_folder,
		"observables": observables, "load_file": load_file,
		"save_to_binary": save_to_binary, "base_parameters": base_parameters,
		"flow_epsilon": flow_epsilon, "NFlows": NFlows,
		"create_perflow_data": create_perflow_data,
		"correct_energy": correct_energy,
		"topct_euclidean_indexes": t_euclidean_indexes,
		"qzero_flow_times": qzero_flow_times}

	databeta60 = copy.deepcopy(default_params)
	databeta60["batch_name"] = beta_folder_name[0]
	databeta60["NCfgs"] = 1000
	databeta60["obs_file"] = "24_6.00"

	databeta61 = copy.deepcopy(default_params)
	databeta61["batch_name"] = beta_folder_name[1]
	databeta61["NCfgs"] = 500
	databeta61["obs_file"] = "28_6.10"

	databeta62 = copy.deepcopy(default_params)
	databeta62["batch_name"] = beta_folder_name[2]
	databeta62["NCfgs"] = 500
	databeta62["obs_file"] = "32_6.20"

	#### Adding relevant batches to args
	# analysis_parameter_list = [databeta60, databeta61, databeta62]
	analysis_parameter_list = [databeta62]

	#### Submitting observable-batches
	for analysis_parameters in analysis_parameter_list:
		analyse(analysis_parameters)

if __name__ == '__main__':
	main()