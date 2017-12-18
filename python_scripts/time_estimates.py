
# TIME ESTIMATES FLOWING
elapsed_time_b60 = (1*24*60*60 + 13*60*60 + 14*60)
NFlows_B60 = 1000
NFlows_B61 = 500
NFlows_B62 = 500
NFlows_B645 = 250

flow_cfg_time_b60 = float(elapsed_time_b60) / float(NFlows_B60)


lattice_size_b60 = 24**3*48.
lattice_size_b61 = 28**3*56.
lattice_size_b62 = 32**3*64.
lattice_size_b645 = 48**3*96.

scaling_b60_b61 = lattice_size_b61/lattice_size_b60
scaling_b60_b62 = lattice_size_b62/lattice_size_b60
scaling_b60_b645 = lattice_size_b645/lattice_size_b60

print "BETA=6.0:\n    %d seconds/%.2f minutes/%.2f hours\n    Time per flow: %.4f minutes" % (elapsed_time_b60,elapsed_time_b60/60.,elapsed_time_b60/3600.,flow_cfg_time_b60/60.0)
print "BETA=6.1:\n    %d hours\n    Time per flow: %.4f minutes" % (elapsed_time_b60/3600. * scaling_b60_b61 * 0.5,flow_cfg_time_b60/60.0 * scaling_b60_b61)
print "BETA=6.2:\n    %d hours\n    Time per flow: %.4f minutes" % (elapsed_time_b60/3600. * scaling_b60_b62 * 0.5,flow_cfg_time_b60/60.0 * scaling_b60_b62)
print "BETA=6.45:\n    %d hours\n    Time per flow: %.4f minutes" % (elapsed_time_b60/3600. * scaling_b60_b645 * 0.25,flow_cfg_time_b60/60.0 * scaling_b60_b645)

import os
# bf = "/work/users/hmvege/"
# ipf = "/work/users/hmvege/output/prodRunBeta6_1/field_configurations/"

bf = "/work/users/hmvege/"
ipf = "//work/users/hmvege/output/prodRunBeta6_1/field_configurations"
ipf = ipf.replace("//","/")

print os.path.normpath(ipf)
print os.path.relpath(ipf,bf)