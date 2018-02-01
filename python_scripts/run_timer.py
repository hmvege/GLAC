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

		# CONFIGURATION GENERATION TIMES [min]
		beta60_cfg = 0.71
		beta61_cfg = 1.48
		beta62_cfg = 2.18

		cfg_times = { 	str(self.beta_values[0]) : beta60_cfg,
						str(self.beta_values[1]) : beta61_cfg,
						str(self.beta_values[2]) : beta62_cfg}

		# Estimates time for 6.45
		b645_cfg_estimates = [(cfg_times[str(beta)]/float(b_corr_updates[str(beta)]))*b_corr_updates[str("6.45")]*b645_size_scaling[str(beta)] for beta in self.beta_values[:3]]
		b645_cfg_estimate = np.mean(b645_cfg_estimates)
		self.b645_cfg_estimate_std = np.std(b645_cfg_estimates)
		cfg_times[str(self.beta_values[3])] = b645_cfg_estimate

		self.cfg_times_per_config = {str(key):cfg_times[str(key)]/float(b_corr_updates[str(key)]) for key in self.beta_values}
		# print {key:val*60 for (key,val) in zip(cfg_times_per_config.keys(),cfg_times_per_config.values())}, "seconds"

		# CONFIGURATION FLOW TIMES [min]
		beta60_flow = 2.3
		beta61_flow = 4.1
		beta62_flow = 8.4

		self.flow_times = {	str(self.beta_values[0]) : beta60_flow,
							str(self.beta_values[1]) : beta61_flow,
							str(self.beta_values[2]) : beta62_flow}

		# Estimates time for 6.45
		b645_flow_estimates = [self.flow_times[str(beta)]*b645_size_scaling[str(beta)] for beta in self.beta_values[:3]]
		b645_flow_estimate = np.mean(b645_flow_estimates)
		self.b645_flow_estimate_std = np.std(b645_flow_estimates)
		self.flow_times[str(self.beta_values[3])] = b645_flow_estimate

	def get_time(self,beta,N_configs,N_corr,N_flows):
		config_time_est = self.cfg_times_per_config[str(beta)]
		total_config_time_est = config_time_est*N_configs*N_corr
		flow_time_est = self.flow_times[str(beta)]
		total_flow_time_est = flow_time_est*N_flows

		error_cfg = ""
		error_flow = ""
		if beta == 6.45:
			error_cfg = "+/- %-.4f hours" % (self.b645_cfg_estimate_std/60.)
			error_flow = "+/- %-.4f hours" % (self.b645_flow_estimate_std/60.)

		print """
Time estimate for run:
beta                %.2f
N_configs           %d
N_corr              %d
N_flows             %d
%s
Time per config:    %10.4f minutes
Total config time:  %10.4f minutes / %-.1f hours %s
Time per flow:      %10.4f minutes
Total flow-time:    %10.4f minutes / %-.1f hours %s
%s""" % (beta,N_configs,N_corr,N_flows,100*"=",
			config_time_est,
			total_config_time_est,total_config_time_est/60,error_flow,
			flow_time_est,
			total_flow_time_est,total_flow_time_est/60,error_flow,
			100*"=")

if __name__ == '__main__':
	description_string = """
Small program intended for estimating times
"""
	parser = argparse.ArgumentParser(prog='Run time estimator', description=description_string)

	######## Program basics #########
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	
	######## Program options ########
	parser.add_argument('beta', 					type=float, help='Beta value to estimate time for.')
	parser.add_argument('-NCfgs','--NConfigs',		type=int, default=None,help='Number of configurations')
	parser.add_argument('-NF','--Nflows',			type=int, default=None,help='Number of flows')
	parser.add_argument('-NCorr','--NCorrelations',	type=int, default=None,help='Number of correlation updates')

	# Default values
	defaults = {"6.0": {"NCorr":200,"NFlows":1000,"NCfgs":1000},
				"6.1": {"NCorr":200,"NFlows":1000,"NCfgs":500},
				"6.2": {"NCorr":200,"NFlows":1000,"NCfgs":500},
				"6.45": {"NCorr":800,"NFlows":1000,"NCfgs":250}}

	# Parses arguments
	if len(sys.argv) == 1:
		args = parser.parse_args(["6.45","-NCfgs","306"])
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

	t = TimeEstimator()
	t.get_time(args.beta,NCfgs,N_corr=NCorr,N_flows=NFlows)
