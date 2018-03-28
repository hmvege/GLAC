import matplotlib.pyplot as plt
import numpy as np
import os

# fig = plt.figure()

# ax1 = fig.add_subplot(221)
# ax2 = fig.add_subplot(222)
# ax3 = fig.add_subplot(223)
# ax4 = fig.add_subplot(224)

# y1 = lambda _x: np.cos(_x**2) + np.random.random(100)
# y2 = lambda _x: np.cosh(np.log(_x**2)) + np.random.random(100)
# y3 = lambda _x: np.sin(_x**2) + np.random.random(100)
# y4 = lambda _x: np.sinh(np.log(_x**2)) + np.random.random(100)

# x = np.linspace(0,10,100)

# ax1.plot(x,y1(x))
# ax2.plot(x,y2(x))
# ax3.plot(x,y3(x))
# ax4.plot(x,y4(x))

# plt.suptitle("Super title!")
# plt.xlabel(r"$x$")
# plt.ylabel(r"$f_{axes}(x)$")

# plt.show()

# topc = np.loadtxt("beta60_correct_topc_flow_config00000.dat", skiprows=3)
# topct = np.loadtxt("beta60_correct_topct_flow_config00000.dat", skiprows=3)

print "Available files:"
print [f for f in os.listdir(os.getcwd()) if os.path.splitext(f)[-1] == ".dat"], "\n"

cfg_num = 95

topc = np.loadtxt("beta60_correct_topc_flow_config00%03d.dat" % cfg_num, skiprows=3)
topct = np.loadtxt("beta60_correct_topct_flow_config00%03d.dat" % cfg_num, skiprows=3)

summed_topc = np.sum(topct[:,1:], axis=1)

eps = 1e-14
diff_checker = lambda _diff: True if _diff < eps else False

for i_topc, i_topct in zip(topc[:,1], summed_topc):
	diff = i_topc - i_topct
	if not diff_checker:
		print "%20.16f %20.16f %20.16f %s" % (i_topc, i_topct, diff, diff_checker(diff)) 
else:
	print "All good: no discrepancies larger than %g" % eps
