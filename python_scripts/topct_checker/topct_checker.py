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

cfg_num = 999
beta = 1
topc = load_file("../../data5/beta6{0:<1d}/flow_observables/topc/prodRunBeta6_{0:<1d}_topc_flow_config00{1:<03d}.dat".format(beta, cfg_num))
topct = load_file("../../data5/beta6{0:<1d}/flow_observables/topct/prodRunBeta6_{0:<1d}_topct_flow_config00{1:<03d}.dat".format(beta, cfg_num))

summed_topc = np.sum(topct[:,1:], axis=1)

eps = 1e-14
diff_checker = lambda _diff: True if np.abs(_diff) > eps else False

# for i_topc, i_topct in zip(topc[:,1], summed_topc):
# 	diff = i_topc - i_topct
# 	if diff_checker(diff):
# 		print "\nDiscrepancies found:"
# 		print "%20.16f %20.16f %20.16f %s" % (i_topc, i_topct, np.abs(diff), diff_checker(diff)) 
# 		break
# else:
# 	print "\nAll good: no discrepancies larger than %g" % eps

diff = np.max(topc[:,1] - summed_topc)
diff_argmax = np.argmax(topc[:,1] - summed_topc)
if diff_checker(diff):
	print "\nDiscrepancies found!"
	print "Max discrepancies:"
	print "%20.16f %20.16f %20.16f %s" % (topc[diff_argmax,1], summed_topc[diff_argmax], np.abs(diff), diff_checker(diff)) 
	# break
else:
	print "\nAll good: no discrepancies larger than %g" % eps
