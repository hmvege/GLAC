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

topc = load_file("beta60_correct_topc_flow_config00%03d.dat" % cfg_num)
topct = load_file("beta60_correct_topct_flow_config00%03d.dat" % cfg_num)

summed_topc = np.sum(topct[:,1:], axis=1)

eps = 1e-14
diff_checker = lambda _diff: True if _diff < eps else False

for i_topc, i_topct in zip(topc[:,1], summed_topc):
	diff = i_topc - i_topct
	if not diff_checker:
		print "\n%20.16f %20.16f %20.16f %s" % (i_topc, i_topct, diff, diff_checker(diff)) 
else:
	print "\nAll good: no discrepancies larger than %g" % eps
