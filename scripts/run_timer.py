import argparse, numpy as np, sys

# All estimates are for 200 correlation updates, so in essence 1 mc update is time / NCorr

class TimeEstimator:
	beta_values = [6.0,6.1,6.2,6.45]

	lattice_size_b60 = 24**3*48.
	lattice_size_b61 = 28**3*56.
	lattice_size_b62 = 32**3*64.
	lattice_size_b645 = 48**3*96.

	NUp10_8x16_time_per_update = 0.082336 / 10 # [seconds / size]
	NUp30_8x16_time_per_update = 0.134728 / 30 # [seconds / size]

	NUScaling = NUp30_8x16_time_per_update / NUp10_8x16_time_per_update

	def __init__(self):
		b645_size_scaling = {	str(self.beta_values[0]) : self.lattice_size_b645 / self.lattice_size_b60,
								str(self.beta_values[1]) : self.lattice_size_b645 / self.lattice_size_b61,
								str(self.beta_values[2]) : self.lattice_size_b645 / self.lattice_size_b62}

		# Correlation updates
		b_corr_updates = {	str(self.beta_values[0]) : 200,
							str(self.beta_values[1]) : 200,
							str(self.beta_values[2]) : 200,
							str(self.beta_values[3]) : 800}

		#### CONFIGURATION GENERATION TIMES [min]
		beta60_cfg = 0.71
		beta61_cfg = 1.48
		beta62_cfg = 2.18

		cfg_times = { 	str(self.beta_values[0]) : beta60_cfg,
						str(self.beta_values[1]) : beta61_cfg,
						str(self.beta_values[2]) : beta62_cfg}

		# Estimates time for 6.45
		b645_cfg_estimates = [(cfg_times[str(beta)]/float(b_corr_updates[str(beta)]))*b_corr_updates[str("6.45")]*b645_size_scaling[str(beta)] for beta in self.beta_values[:3]]
		self.b645_cfg_estimate_std = np.std(b645_cfg_estimates) / b_corr_updates[str("6.45")]
		cfg_times[str(self.beta_values[3])] = np.mean(b645_cfg_estimates)

		# Normalizes by N updates
		self.cfg_times_per_config = {str(key):cfg_times[str(key)]/float(b_corr_updates[str(key)]) for key in self.beta_values}
		# print {key:val*60 for (key,val) in zip(cfg_times_per_config.keys(),cfg_times_per_config.values())}, "seconds"

		#### CONFIGURATION FLOW TIMES [min]
		beta60_flow = 2.3
		beta61_flow = 4.1
		beta62_flow = 8.4

		# Time per flow per config
		self.flow_times = {	str(self.beta_values[0]) : beta60_flow / 1000,
							str(self.beta_values[1]) : beta61_flow / 1000,
							str(self.beta_values[2]) : beta62_flow / 1000}

		# Estimates time for 6.45
		b645_flow_estimates = [self.flow_times[str(beta)]*b645_size_scaling[str(beta)] for beta in self.beta_values[:3]]

		self.b645_flow_estimate_std = np.std(b645_flow_estimates)
		self.flow_times[str(self.beta_values[3])] = np.mean(b645_flow_estimates)

	def get_time(self, beta, N_configs, N_corr, N_therms, N_flows, numprocs=512, N_updates=10, verbose=False):
		NUpdatesScaling = 1 + self.NUScaling * (N_updates - 10) / 20

		# Test runs done with 512 cores, scaling around 0.55
		cpu_scaling = (0.55 ** ((numprocs-512) / 512))

		config_time_est = self.cfg_times_per_config[str(beta)] * NUpdatesScaling

		total_config_time_est = config_time_est * N_configs * N_corr 
		config_cpu_time = total_config_time_est*cpu_scaling/60.0*numprocs

		# Estimates thermalization
		thermalization_time_est = config_time_est * N_therms

		# Estimates flow times
		flow_time_est = self.flow_times[str(beta)]
		total_flow_time_est = flow_time_est*N_configs*N_flows
		flow_cpu_time = total_flow_time_est*cpu_scaling/60.0*numprocs

		# Total time for program run
		total_time = total_flow_time_est + total_config_time_est + thermalization_time_est


		total_cpu_time = total_time*cpu_scaling/60.0*numprocs

		if verbose:
			print """
N updates scaling factor:   %d

""" % (NUpdatesScaling)

		error_cfg = ""
		error_therm = ""
		error_flow = ""
		total_error = ""
		if beta == 6.45:
			num_error_cfg = self.b645_cfg_estimate_std/60.*N_configs*N_corr
			num_error_flow = self.b645_flow_estimate_std/60.*N_configs*N_flows
			num_error_therm = self.b645_cfg_estimate_std/60.*N_therms
			error_cfg = "+/- %-.6f hours" % (num_error_cfg)
			error_therm = "+/- %-.6f hours" % (num_error_therm)
			error_flow = "+/- %-.6f hours" % (num_error_flow)
			total_error = "+/- %-.6f hours" % (np.sqrt(num_error_cfg**2 + num_error_flow**2))

		msg = "Time estimate for run:"
		msg += "\nbeta                %.2f" % beta
		msg += "\nN_configs           %d" % N_configs
		msg += "\nN_corr              %d" % N_corr
		msg += "\nN_therm             %d" % N_therms
		msg += "\nN_updates           %d" % N_updates
		msg += "\nN_flows             %d" % N_flows
		msg += "\nCPUs                %d" % numprocs
		msg += "\n" + 100*"="
		msg += "\nTime per config:    %10.3f minutes" % (config_time_est*cpu_scaling)
		msg += "\nTotal config time:  %10.2f minutes / %-.1f hours %s" % (total_config_time_est*cpu_scaling, total_config_time_est*cpu_scaling/60, error_cfg)
		msg += "\nConfig CPU hours:   %10.2f hours" % config_cpu_time
		msg += "\nThermalization time:%10.2f minutes / %-.1f hours %s" % (thermalization_time_est*cpu_scaling, thermalization_time_est*cpu_scaling/60, error_therm)
		msg += "\nTime per flow:      %10.3f minutes" % (flow_time_est*cpu_scaling)
		msg += "\nTotal flow-time:    %10.2f minutes / %-.1f hours %s" % (total_flow_time_est*cpu_scaling, total_flow_time_est*cpu_scaling/60, error_flow)
		msg += "\nFlow CPU hours:     %10.2f hours" % flow_cpu_time
		msg += "\nTotal time:         %10.2f minutes / %-.1f hours %s" % (total_time*cpu_scaling, (total_time*cpu_scaling)/60, total_error)
		msg += "\nTotal CPU hours:    %10.2f hours" % total_cpu_time
		msg += "\n" + 100*"="

		print msg

if __name__ == '__main__':
	description_string = """Small program intended for estimating times"""

	parser = argparse.ArgumentParser(prog='Run time estimator', description=description_string)

	######## Program basics #########
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	
	######## Program options ########
	parser.add_argument('beta', 						type=float, help='Beta value to estimate time for.')
	parser.add_argument('-NCfgs','--NConfigs',			type=int, default=None,help='Number of configurations')
	parser.add_argument('-NF','--Nflows',				type=int, default=None,help='Number of flows')
	parser.add_argument('-NCorr','--NCorrelations',		type=int, default=None,help='Number of correlation updates')
	parser.add_argument('-NTherm','--NThermalizations',	type=int, default=None,help='Number of thermalization updates')
	parser.add_argument('-NUp','--NUpdates',			type=int, default=None,help='Number of updates per link')
	parser.add_argument('-numprocs',					type=int, default=None,help='Number of processors')

	# Default values
	defaults = {"6.0": {"NCorr":200,"NFlows":1000,"NCfgs":1000,"NTherm":20000,"NUpdates":10,"numprocs":512},
				"6.1": {"NCorr":200,"NFlows":1000,"NCfgs":500,"NTherm":20000,"NUpdates":10,"numprocs":512},
				"6.2": {"NCorr":200,"NFlows":1000,"NCfgs":500,"NTherm":20000,"NUpdates":10,"numprocs":512},
				"6.45": {"NCorr":800,"NFlows":1000,"NCfgs":250,"NTherm":20000,"NUpdates":10,"numprocs":512}}

	# Parses arguments
	# if len(sys.argv) == 1:
	# args = parser.parse_args(['6.1', '-NCorr', '600', '-NUp', '30', '-NF', '1000', '-numprocs', '512'])
	# 	args = parser.parse_args(["6.45","-NCfgs","250"])
	args = parser.parse_args()
	# args = parser.parse_args(['6.0', '-NCorr', '600', '-NCfgs', '1000', '-NF', '1000', '-NTherm', '20000', '-NUp', '30'])
	if str(args.beta) not in defaults.keys():
		raise KeyError("Error: valid beta values: %s" % ", ".join(sorted(defaults.keys())))

	# Checks and retrives configurations
	if args.NConfigs != None:
		NCfgs = args.NConfigs
	else:
		NCfgs = defaults[str(args.beta)]["NCfgs"]

	# Checks and retrives flows
	if args.Nflows != None:
		NFlows = args.Nflows
	else:
		NFlows = defaults[str(args.beta)]["NFlows"]

	# Checks and retrives correlation updates
	if args.NCorrelations != None:
		NCorr = args.NCorrelations
	else:
		NCorr = defaults[str(args.beta)]["NCorr"]

	# Checks and retrieves thermalization updates
	if args.NThermalizations != None:
		NTherm = args.NThermalizations
	else:
		NTherm = defaults[str(args.beta)]["NTherm"]		

	if args.NUpdates != None:
		NUpdates = args.NUpdates
	else:
		NUpdates = defaults[str(args.beta)]["NUpdates"]

	if args.numprocs != None:
		numprocs = args.numprocs
	else:
		numprocs = defaults[str(args.beta)]["numprocs"]

	t = TimeEstimator()
	t.get_time(args.beta,NCfgs,N_corr=NCorr,N_therms= NTherm,N_flows=NFlows,numprocs=numprocs,N_updates=NUpdates)
