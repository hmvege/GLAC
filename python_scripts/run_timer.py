import argparse, numpy as np, sys

# All estimates are for 200 correlation updates, so in essence 1 mc update is time / NCorr

class TimeEstimator:
	beta_values = [6.0,6.1,6.2,6.45]

	lattice_size_b60 = 24**3*48.
	lattice_size_b61 = 28**3*56.
	lattice_size_b62 = 32**3*64.
	lattice_size_b645 = 48**3*96.

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

	def get_time(self,beta,N_configs,N_corr,N_therms,N_flows,numprocs=512):
		# Estimates configuration times
		config_time_est = self.cfg_times_per_config[str(beta)]
		total_config_time_est = config_time_est*N_configs*N_corr

		# Estimates thermalization
		thermalization_time_est = config_time_est * N_therms

		# Estimates flow times
		flow_time_est = self.flow_times[str(beta)]
		total_flow_time_est = flow_time_est*N_configs*N_flows

		total_time = total_flow_time_est + total_config_time_est + thermalization_time_est

		# Test runs done with 512 cores, scaling around 0.55
		cpu_scaling = (0.55 ** ((numprocs-512) / 512))
		print cpu_scaling
		cpu_time = total_time*cpu_scaling/60.0*numprocs

		error_cfg = ""
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

		print """
Time estimate for run:
beta                %.2f
N_configs           %d
N_corr              %d
N_therms            %d
N_flows             %d
CPUs                %d
%s
Time per config:    %10.3f minutes
Total config time:  %10.2f minutes / %-.1f hours %s
Thermalization time:%10.2f minutes / %-.1f hours %s
Time per flow:      %10.3f minutes
Total flow-time:    %10.2f minutes / %-.1f hours %s
Total time:         %10.2f minutes / %-.1f hours %s
CPU hours:          %10.2f hours
%s""" % (beta,N_configs,N_corr,N_therms,N_flows,numprocs,100*"=",
			config_time_est*cpu_scaling,
			total_config_time_est*cpu_scaling,total_config_time_est*cpu_scaling/60,error_cfg,
			thermalization_time_est*cpu_scaling,thermalization_time_est*cpu_scaling/60,error_therm,
			flow_time_est*cpu_scaling,
			total_flow_time_est*cpu_scaling,total_flow_time_est*cpu_scaling/60,error_flow,
			total_time*cpu_scaling, (total_time*cpu_scaling)/60, total_error,
			cpu_time,
			100*"=")

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
	parser.add_argument('-numprocs',					type=int, default=None,help='Number of processors')

	# Default values
	defaults = {"6.0": {"NCorr":200,"NFlows":1000,"NCfgs":1000,"NTherm":20000,"numprocs":512},
				"6.1": {"NCorr":200,"NFlows":1000,"NCfgs":500,"NTherm":20000,"numprocs":512},
				"6.2": {"NCorr":200,"NFlows":1000,"NCfgs":500,"NTherm":20000,"numprocs":512},
				"6.45": {"NCorr":800,"NFlows":1000,"NCfgs":250,"NTherm":20000,"numprocs":512}}

	# Parses arguments
	if len(sys.argv) == 1:
		args = parser.parse_args(["6.45","-NCfgs","250"])
	else:
		args = parser.parse_args()

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

	if args.numprocs != None:
		numprocs = args.numprocs
	else:
		numprocs = defaults[str(args.beta)]["numprocs"]

	t = TimeEstimator()
	t.get_time(args.beta,NCfgs,N_corr=NCorr,N_therms= NTherm,N_flows=NFlows,numprocs=numprocs)
