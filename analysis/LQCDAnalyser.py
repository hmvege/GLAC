from pre_analysis.flow_analyser import *
from tools.folderreadingtools import DataReader
from tools.postanalysisdatareader import PostAnalysisDataReader
from post_analysis.post_analysis import *
import statistics.parallel_tools as ptools
import os
import numpy as np
import copy
import sys
import time

def analyse_default(analysis_object, N_bs, NBins=30):
	print analysis_object
	analysis_object.boot(N_bs)
	analysis_object.jackknife()
	analysis_object.save_post_analysis_data()
	analysis_object.plot_original()
	analysis_object.plot_boot()
	analysis_object.plot_jackknife()
	analysis_object.autocorrelation()
	analysis_object.plot_autocorrelation(0)
	analysis_object.plot_autocorrelation(-1)
	analysis_object.plot_mc_history(0)
	analysis_object.plot_mc_history(int(analysis_object.NFlows * 0.25))
	analysis_object.plot_mc_history(int(analysis_object.NFlows * 0.50))
	analysis_object.plot_mc_history(int(analysis_object.NFlows * 0.75))
	analysis_object.plot_mc_history(-1)
	analysis_object.plot_original()
	analysis_object.plot_boot()
	analysis_object.plot_jackknife()
	analysis_object.plot_histogram(0, NBins=NBins)
	analysis_object.plot_histogram(int(analysis_object.NFlows * 0.25), NBins=NBins)
	analysis_object.plot_histogram(int(analysis_object.NFlows * 0.50), NBins=NBins)
	analysis_object.plot_histogram(int(analysis_object.NFlows * 0.75), NBins=NBins)
	analysis_object.plot_histogram(-1, NBins=NBins)
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

	analyse_default(topc_analysis, N_bs, NBins=150)

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

def analyse_qtq0e(params, flow_time_indexes, euclidean_time_percents):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	qtqzero_analysis = AnalyseQtQ0Euclidean(obs_data("topct"), dryrun=dryrun,
		parallel=parallel, numprocs=numprocs, verbose=verbose)

	for flow_time_index in flow_time_indexes:
		for euclidean_percent in euclidean_time_percents:
			qtqzero_analysis.set_flow_time(flow_time_index, euclidean_percent)
			analyse_default(qtqzero_analysis, N_bs)

def analyse_topct(params, numsplits):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	topct_analysis = AnalyseTopChargeInEuclTime(obs_data("topct"),
		dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)

	indexes = np.linspace(0, topct_analysis.NT, numsplits, dtype=int) - 1
	indexes[0] += 1
	for ie in indexes:
		topct_analysis.setEQ0(ie)
		analyse_default(topct_analysis, N_bs)

def analyse_topsust(params, numsplits):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	topct_analysis = AnalyseTopSusInEuclTime(obs_data("topct"),
		dryrun=dryrun, parallel=parallel, numprocs=numprocs, verbose=verbose)

	indexes = np.linspace(0, topct_analysis.NT, numsplits, dtype=int) - 1
	indexes[0] += 1
	for ie in indexes:
		topct_analysis.setEQ0(ie)
		analyse_default(topct_analysis, N_bs)

def analyse_topcEuclTime(params, numsplits=None, intervals=None):
	"""
	Analysis function for the topological charge in euclidean time. Requires
	either numsplits or intervals.

	Args:
		params: list of default parameters containing obs_data, dryrun, 
			parallel, numprocs, verbose, N_bs.
		numsplits: number of splits to make in the dataset. Default is None.
		intervals: intervals to plot in. Default is none.
	"""
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	analyse_topcte = AnalyseTopChargeSplitEuclideanTime(obs_data("topct"),
		dryrun=dryrun, parallel=parallel, 
		numprocs=numprocs, verbose=verbose)

	# Sets up the intervals
	if (intervals == numsplits == None) or (intervals != None and numsplits != None):
		raise KeyError("Either provide intervals to plot for or the number of intervals to split into.")

	NT = analyse_topcte.NT
	if intervals == None:
		split_interval = NT/numsplits
		intervals = zip(
			range(0, NT+1, split_interval), 
			range(split_interval, NT+1, split_interval)
		)
		assert NT % numsplits == 0, "Bad number of splits: NT % numplits = %d " % (NT % numsplits)

	t_interval = iter(intervals)

	for t_int in t_interval:
		analyse_topcte.set_t_interval(t_int)
		analyse_default(analyse_topcte, N_bs)

def analyse_topsusEuclTime(params, numsplits=None, intervals=None):
	"""
	Analysis function for the topological susceptibility in euclidean time. 
	Requires either numsplits or intervals.

	Args:
		params: list of default parameters containing obs_data, dryrun, 
			parallel, numprocs, verbose, N_bs.
		numsplits: number of splits to make in the dataset. Default is None.
		intervals: intervals to plot in. Default is none.
	"""
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	analyse_topsuste = AnalysisTopSusSplitEuclideanTime(obs_data("topct"),
		dryrun=dryrun, parallel=parallel, 
		numprocs=numprocs, verbose=verbose)

	# Sets up the intervals
	if (intervals == numsplits == None) or (intervals != None and numsplits != None):
		raise KeyError("Either provide intervals to plot for or the number of intervals to split into.")

	NT = analyse_topsuste.NT
	if intervals == None:
		split_interval = NT/numsplits
		intervals = zip(
			range(0, NT+1, split_interval), 
			range(split_interval, NT+1, split_interval)
		)
		assert NT % numsplits == 0, "Bad number of splits: NT % numplits = %d " % (NT % numsplits)

	t_interval = iter(intervals)

	for t_int in t_interval:
		analyse_topsuste.set_t_interval(t_int)
		analyse_default(analyse_topsuste, N_bs)


def analyse_topcMCTime(params, numsplits=None, intervals=None):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	analyse_topcMC = AnalysisTopChargeSplitMCTime(obs_data("topc"),
		dryrun=dryrun, parallel=parallel,
		numprocs=numprocs, verbose=verbose)

	# Sets up the intervals
	if (intervals == numsplits == None) or (intervals != None and numsplits != None):
		raise KeyError("Either provide MC intervals to plot for or the number of MC intervals to split into.")

	NCfgs = analyse_topcMC.N_configurations
	if intervals == None:
		split_interval = NCfgs/numsplits
		intervals = zip(
			range(0, NCfgs+1, split_interval), 
			range(split_interval, NCfgs+1, split_interval)
		)
		assert NCfgs % numsplits == 0, "Bad number of splits: NCfgs % numplits = %d " % (NCfgs % numsplits)

	MC_interval = iter(intervals)

	for MC_int in MC_interval:
		analyse_topcMC.set_MC_interval(MC_int)
		analyse_default(analyse_topcMC, N_bs)

def analyse_topsusMCTime(params, numsplits=None, intervals=None):
	obs_data, dryrun, parallel, numprocs, verbose, N_bs = params

	analyse_topcMC = AnalysisTopSusSplitMCTime(obs_data("topc"),
		dryrun=dryrun, parallel=parallel,
		numprocs=numprocs, verbose=verbose)

	# Sets up the intervals
	if (intervals == numsplits == None) or (intervals != None and numsplits != None):
		raise KeyError("Either provide MC intervals to plot for or the number of MC intervals to split into.")

	NCfgs = analyse_topcMC.N_configurations
	if intervals == None:
		split_interval = NCfgs/numsplits
		intervals = zip(
			range(0, NCfgs+1, split_interval), 
			range(split_interval, NCfgs+1, split_interval)
		)
		assert NCfgs % numsplits == 0, "Bad number of splits: NCfgs % numplits = %d " % (NCfgs % numsplits)

	MC_interval = iter(intervals)

	for MC_int in MC_interval:
		analyse_topcMC.set_MC_interval(MC_int)
		analyse_default(analyse_topcMC, N_bs)


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
	print "="*100

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
	if "qtq0e" in parameters["observables"]:
		analyse_qtq0e(params, parameters["flow_time_indexes"], parameters["euclidean_time_percents"])
	if "topc4" in parameters["observables"]:
		analyse_topc4(params)
	if "topsust" in parameters["observables"]:
		analyse_topsust(params, parameters["num_t_euclidean_indexes"])
	if "topct" in parameters["observables"]:
		analyse_topct(params, parameters["num_t_euclidean_indexes"])
	if "topcte" in parameters["observables"]:
		analyse_topcEuclTime(params, parameters["numsplits_eucl"], parameters["intervals_eucl"])
	if "topsuste" in parameters["observables"]:
		analyse_topsusEuclTime(params, parameters["numsplits_eucl"], parameters["intervals_eucl"])
	if "topcMC" in parameters["observables"]:
		analyse_topcMCTime(params, parameters["MC_time_splits"])
	if "topsusMC" in parameters["observables"]:
		analyse_topsusMCTime(params, parameters["MC_time_splits"])
	
	post_time = time.clock()
	print "="*100
	print "Analysis of batch %s observables %s in %.2f seconds" % (batch_name,
		", ".join([i.lower() for i in parameters["observables"]]), (post_time-pre_time))
	print "="*100

def post_analysis(batch_folder, batch_beta_names, observables, topsus_fit_target,
	line_fit_interval, energy_fit_target, post_analysis_data_type=None, 
	bval_to_plot="all", verbose=False):
	"""
	Post analysis of the flow observables.

	Args: 
		batch_folder: string, folder containing all the beta data.
		batch_beta_names: list of the beta folder names.
		topsus_fit_target: list of x-axis points to line fit at.
		line_fit_interval: float, extension of the area around the fit target 
			that will be used for the line fit.
		energy_fit_target: point of which we will perform a line fit at.
	"""

	print "="*100 + "\nPost-analysis: retrieving data from: %s" % batch_folder

	if post_analysis_data_type == None:
		post_analysis_data_type = ["bootstrap", "jackknife"]

	# Loads data from post analysis folder
	data = PostAnalysisDataReader(batch_folder)

	for analysis_type in post_analysis_data_type:
		if "topsus" in observables:
			# Plots topsus
			topsus_analysis = TopSusPostAnalysis(data, verbose=verbose)
			print topsus_analysis
			topsus_analysis.set_analysis_data_type(analysis_type)
			topsus_analysis.plot()

			# Retrofits the topsus for continuum limit
			continium_targets = [0.3, 0.4, 0.5, 0.58]
			for cont_target in continium_targets:
				topsus_analysis.plot_continuum(cont_target)

		if "energy" in observables:
			# Plots energy
			energy_analysis = EnergyPostAnalysis(data, verbose=verbose)
			print energy_analysis
			energy_analysis.set_analysis_data_type(analysis_type)
			energy_analysis.plot()

			# # Retrofits the energy for continiuum limit
			# energy_analysis.plot_continuum(0.3, 0.015, "bootstrap_fit")

			# # Plot running coupling
			# energy_analysis.coupling_fit()

		if "topc" in observables:
			topc_analysis = TopChargePostAnalysis(data, verbose=verbose)
			print topc_analysis
			topc_analysis.set_analysis_data_type(analysis_type)
			topc_analysis.plot(y_limits=[-5,5])

		if "plaq" in observables:
			plaq_analysis = PlaqPostAnalysis(data, verbose=verbose)
			print plaq_analysis
			plaq_analysis.set_analysis_data_type(analysis_type)
			plaq_analysis.plot()

		if "topc4" in observables:
			topc4_analysis = TopCharge4PostAnalysis(data, verbose=verbose)
			print topc4_analysis
			topc4_analysis.set_analysis_data_type(analysis_type)
			topc4_analysis.plot()

		if "qtqzero" in observables:
			qtq0_analysis = QtQ0PostAnalysis(data, verbose=verbose)
			print qtq0_analysis
			qtq0_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = qtq0_analysis.get_N_intervals()
			for i in range(N_int):
				qtq0_analysis.plot_interval(i)
			qtq0_analysis.plot_series([0,1,2,3], beta=bval_to_plot)
			qtq0_analysis.plot_series([4,5,6,7], beta=bval_to_plot)

		if "topct" in observables:
			topct_analysis = TopChargeTPostAnalysis(data, verbose=verbose)
			print topct_analysis
			topct_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = topct_analysis.get_N_intervals()
			for i in range(N_int):
				topct_analysis.plot_interval(i)
			topct_analysis.plot_series([0,1,2,3], beta=bval_to_plot)

		if "topcte" in observables:
			topcte_analysis = TopChargeEuclSplitPostAnalysis(data, verbose=verbose)
			print topcte_analysis
			topcte_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = topcte_analysis.get_N_intervals()
			for i in range(N_int):
				topcte_analysis.plot_interval(i)
			topcte_analysis.plot_series([0,1,2,3], beta=bval_to_plot)

		if "topcMC" in observables:
			topcmc_analysis = TopChargeMCSplitPostAnalysis(data, verbose=verbose)
			print topcmc_analysis
			topcmc_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = topcmc_analysis.get_N_intervals()
			for i in range(N_int):
				topcmc_analysis.plot_interval(i)
			topcmc_analysis.plot_series([0,1,2,3], beta=bval_to_plot)

		if "topsust" in observables:
			topsust_analysis = TopSusTPostAnalysis(data, verbose=verbose)
			print topsust_analysis
			topsust_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = topsust_analysis.get_N_intervals()
			for i in range(N_int):
				topsust_analysis.plot_interval(i)
			topsust_analysis.plot_series([0,1,2,3], beta=bval_to_plot)

		if "topsuste" in observables:
			topsuste_analysis = TopSusEuclSplitPostAnalysis(data, verbose=verbose)
			print topsuste_analysis
			topsuste_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = topsuste_analysis.get_N_intervals()
			for i in range(N_int):
				topsuste_analysis.plot_interval(i)
			topsuste_analysis.plot_series([0,1,2,3], beta=bval_to_plot)

		if "topsusMC" in observables:
			topsusmc_analysis = TopSusMCSplitPostAnalysis(data, verbose=verbose)
			print topsusmc_analysis
			topsusmc_analysis.set_analysis_data_type(analysis_type)
			N_int, intervals = topsusmc_analysis.get_N_intervals()
			for i in range(N_int):
				topsusmc_analysis.plot_interval(i)
			topsusmc_analysis.plot_series([0,1,2,3], beta=bval_to_plot)


def main():
	#### Available observables
	observables = [
		"plaq", "energy", "topc", "topsus", "qtqzero", "topc4", 
		"topct", "topcte", "topcMC", "topsuste", 
		"topsusMC", "topsust", "qtq0e"
	]

	# ### Basic observables
	# observables = ["plaq", "energy", "topc", "topsus"]
	# observables = basic_observables
	# observables = ["topc", "plaq", "topsus", "topc4"]

	# observables = ["topcMC", "topsusMC"]
	# observables = ["topcte", "topcMC", "topsuste", "topsusMC", "topct"]
	# observables = ["topsuste", "topsusMC", "topct"]
	# observables = ["topsusMC"]
	observables = ["qtq0e"]

	print 100*"=" + "\nObservables to be analysed: %s" % ", ".join(observables)
	print 100*"=" + "\n"

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
	NFlows = 50
	flow_epsilon = 0.01

	#### Post analysis parameters
	run_post_analysis = True
	line_fit_interval = 0.015
	topsus_fit_targets = [0.3,0.4,0.5,0.58]
	energy_fit_target = 0.3

	#### Different batches
	# data_batch_folder = "data2"
	# data_batch_folder = "data4"
	data_batch_folder = "data5"
	# data_batch_folder = "DataGiovanni"
	# data_batch_folder = "smaug_data_beta61"

	#### If we need to multiply
	if data_batch_folder == "DataGiovanni":
		if "topct" in observables:
			observables.remove("topct")
		correct_energy = False
		load_file = True
		save_to_binary = False
	else:
		correct_energy = True

	#### Different beta values folders:
	# beta_folders = ["beta60", "beta61", "beta62"]
	beta_folders = ["beta60", "beta61", "beta62", "beta645"]
	# beta_folders = ["beta6_0", "beta6_1", "beta6_2"]
	# beta_folders = ["beta61"]

	# Indexes to look at for topct.
	num_t_euclidean_indexes = 5

	# Number of different sectors we will analyse in euclidean time
	numsplits_eucl = 4
	intervals_eucl = None

	# Number of different sectors we will analyse in monte carlo time
	MC_time_splits = 4

	# Percents of data where we do qtq0
	qzero_flow_times = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]

	# Flow time indexes to plot qtq0 in euclidean time at
	flow_time_indexes = [0, 50, 200, 400, 700, 999]
	euclidean_time_percents = [0, 0.25, 0.50, 0.75, 1.00]

	#### Analysis batch setups
	default_params = {
		"batch_folder": data_batch_folder,
		"observables": observables, "load_file": load_file,
		"save_to_binary": save_to_binary, "base_parameters": base_parameters,
		"flow_epsilon": flow_epsilon, "NFlows": NFlows,
		"create_perflow_data": create_perflow_data,
		"correct_energy": correct_energy,
		"num_t_euclidean_indexes": num_t_euclidean_indexes,
		"qzero_flow_times": qzero_flow_times,
		"flow_time_indexes": flow_time_indexes,
		"euclidean_time_percents": euclidean_time_percents,
		"numsplits_eucl": numsplits_eucl,
		"intervals_eucl": intervals_eucl,
		"MC_time_splits": MC_time_splits,
	}

	databeta60 = copy.deepcopy(default_params)
	databeta60["batch_name"] = beta_folders[0]
	databeta60["NCfgs"] = 1000
	databeta60["obs_file"] = "24_6.00"

	databeta61 = copy.deepcopy(default_params)
	databeta61["batch_name"] = beta_folders[1]
	databeta61["NCfgs"] = 500
	databeta61["obs_file"] = "28_6.10"

	databeta62 = copy.deepcopy(default_params)
	databeta62["batch_name"] = beta_folders[2]
	databeta62["NCfgs"] = 500
	databeta62["obs_file"] = "32_6.20"

	databeta645 = copy.deepcopy(default_params)
	databeta645["batch_name"] = beta_folders[3]
	databeta645["NCfgs"] = 250
	databeta645["obs_file"] = "48_6.45"


	# smaug_data_beta61_analysis = copy.deepcopy(default_params)
	# smaug_data_beta61_analysis["batch_name"] = beta_folders[0]
	# smaug_data_beta61_analysis["NCfgs"] = 100

	#### Adding relevant batches to args
	analysis_parameter_list = [databeta60, databeta61, databeta62, databeta645]
	analysis_parameter_list = [databeta60, databeta61, databeta62]
	# analysis_parameter_list = [databeta62]
	# analysis_parameter_list = [databeta61, databeta62]
	# analysis_parameter_list = [smaug_data_beta61_analysis]

	#### Submitting observable-batches
	for analysis_parameters in analysis_parameter_list:
		analyse(analysis_parameters)

	# #### Submitting post-analysis data
	# if len(analysis_parameter_list) >= 2:
	# 	post_analysis(data_batch_folder, beta_folders, observables, topsus_fit_targets,
	# 		line_fit_interval, energy_fit_target, verbose=verbose)

if __name__ == '__main__':
	main()