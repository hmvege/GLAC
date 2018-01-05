# TIME ESTAMTES GENERATING
config_elapsed_time_b60 = 42817.5 # seconds
time_per_config_generation = config_elapsed_time_b60 / 1000.0

# TIME ESTIMATES FLOWING
flow_elapsed_time_b60 = (1*24*60*60 + 13*60*60 + 14*60)
NFlows_B60 = 1000
NFlows_B61 = 500
NFlows_B62 = 500
NFlows_B645 = 250

flow_cfg_time_b60 = float(flow_elapsed_time_b60) / float(NFlows_B60)

lattice_size_b60 = 24**3*48.
lattice_size_b61 = 28**3*56.
lattice_size_b62 = 32**3*64.
lattice_size_b645 = 48**3*96.

scaling_b60_b61 = lattice_size_b61/lattice_size_b60
scaling_b60_b62 = lattice_size_b62/lattice_size_b60
scaling_b60_b645 = lattice_size_b645/lattice_size_b60

# print "BETA=6.0:\n    %d seconds/%.2f minutes/%.2f hours\n    Time per flow: %.4f minutes" % (flow_elapsed_time_b60,flow_elapsed_time_b60/60.,flow_elapsed_time_b60/3600.,flow_cfg_time_b60/60.0)
# print "BETA=6.1:\n    %d hours\n    Time per flow: %.4f minutes" % (flow_elapsed_time_b60/3600. * scaling_b60_b61 * 0.5,flow_cfg_time_b60/60.0 * scaling_b60_b61)
# print "BETA=6.2:\n    %d hours\n    Time per flow: %.4f minutes" % (flow_elapsed_time_b60/3600. * scaling_b60_b62 * 0.5,flow_cfg_time_b60/60.0 * scaling_b60_b62)
# print "BETA=6.45:\n    %d hours\n    Time per flow: %.4f minutes" % (flow_elapsed_time_b60/3600. * scaling_b60_b645 * 0.25,flow_cfg_time_b60/60.0 * scaling_b60_b645)

def print_flow_time_estimates(beta, tot_gen_time, time_per_generation, tot_flow_time, time_per_flow_time, scaling=1.0, config_number_scaling_factor=1.0):
	msg = """BETA={0:<g}:
    Total config generation time:   {1:<.0f} seconds/{2:<.1f} minutes/{3:<.2f} hours
    Time per configuration:         {4:<.4f} minutes
    Total flow time: 				{5:<.0f} seconds/{6:<.1f} minutes/{7:<.2f} hours
    Time per flow:   				{8:<.4f} minutes""".format(
		beta,
		tot_gen_time * scaling * config_number_scaling_factor,
		tot_gen_time/60.0 * scaling * config_number_scaling_factor,
		tot_gen_time/3600.0 * scaling * config_number_scaling_factor,
		time_per_generation/60.0 * scaling,
		tot_flow_time * scaling * config_number_scaling_factor,
		tot_flow_time/60.0 * scaling * config_number_scaling_factor,
		tot_flow_time/3600.0 * scaling * config_number_scaling_factor,
		time_per_flow_time/60.0 * scaling)
	print msg

print_flow_time_estimates(6.0, config_elapsed_time_b60, time_per_config_generation, flow_elapsed_time_b60, flow_cfg_time_b60)
print_flow_time_estimates(6.1, config_elapsed_time_b60, time_per_config_generation, flow_elapsed_time_b60, flow_cfg_time_b60,scaling=scaling_b60_b61,config_number_scaling_factor=0.5)
print_flow_time_estimates(6.2, config_elapsed_time_b60, time_per_config_generation, flow_elapsed_time_b60, flow_cfg_time_b60,scaling=scaling_b60_b62,config_number_scaling_factor=0.5)
print_flow_time_estimates(6.45, config_elapsed_time_b60, time_per_config_generation, flow_elapsed_time_b60, flow_cfg_time_b60,scaling=scaling_b60_b645,config_number_scaling_factor=0.25)


# import os
# # bf = "/work/users/hmvege/"
# # ipf = "/work/users/hmvege/output/prodRunBeta6_1/field_configurations/"

# bf = "/work/users/hmvege/"
# ipf = "//work/users/hmvege/output/prodRunBeta6_1/field_configurations"
# ipf = ipf.replace("//","/")

# print os.path.normpath(ipf)
# print os.path.relpath(ipf,bf)