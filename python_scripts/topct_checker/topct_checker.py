import matplotlib.pyplot as plt
import numpy as np
import os

def load_file(fname):
	data = np.loadtxt(fname, skiprows=3)
	print "%s loaded." % fname
	return data

print "Available files:"
print [f for f in os.listdir(os.getcwd()) if os.path.splitext(f)[-1] == ".dat"], "\n"

cfg_num = 352

# topc = load_file("beta60_correct_topc_flow_config00%03d.dat" % cfg_num)
# topct = load_file("beta60_correct_topct_flow_config00%03d.dat" % cfg_num)

cfg_num = 3
beta = 45
eps = 1e-14

# topc_fname = "../../data5/beta6{0:<1d}/flow_observables/topc/prodRunBeta6_{0:<1d}_topc_flow_config00{1:<03d}.dat".format(beta, cfg_num)
# topct_fname = "../../data5/beta6{0:<1d}/flow_observables/topct/prodRunBeta6_{0:<1d}_topct_flow_config00{1:<03d}.dat".format(beta, cfg_num)

## From corrected data8 folder
# topc_fname = "../../data8/beta6{0:<1d}/flow_observables/topc/beta6{0:<1d}_correct_topc_flow_config00{1:>03d}.dat".format(beta, cfg_num)
# topct_fname = "../../data8/beta6{0:<1d}/flow_observables/topct/beta6{0:<1d}_correct_topct_flow_config00{1:>03d}.dat".format(beta, cfg_num)

## From production run data8 folder, mainly beta 6.45
topc_fname = "../../data8/beta6{0:<1d}/flow_observables/topc/prodRunBeta6_{0:<1d}_topc_flow_config00{1:>03d}.dat".format(beta, cfg_num)
topct_fname = "../../data8/beta6{0:<1d}/flow_observables/topct/prodRunBeta6_{0:<1d}_topct_flow_config00{1:>03d}.dat".format(beta, cfg_num)

topc = load_file(topc_fname)
topct = load_file(topct_fname)

print topc.shape
print topct.shape

summed_topc = np.sum(topct[:,1:], axis=1)

diff_checker = lambda _diff: True if np.abs(_diff) > eps else False

diff = np.max(topc[:,1] - summed_topc)
diff_argmax = np.argmax(topc[:,1] - summed_topc)

if diff_checker(diff):
	print "\nDiscrepancies found!"
	print "Max discrepancies:"
	print "%20.16f %20.16f %20.16f %s" % (topc[diff_argmax,1], summed_topc[diff_argmax], np.abs(diff), diff_checker(diff)) 
	# break
else:
	print "\nAll good: no discrepancies larger than %g" % eps
