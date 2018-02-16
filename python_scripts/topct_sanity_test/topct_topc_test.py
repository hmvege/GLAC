import numpy as np, os, matplotlib.pyplot as plt

print os.getcwd(), os.listdir(os.getcwd())

topct_file_path = "../../output/topctTestRun/flow_observables/topct/topctTestRun_topct_flow_config00000.dat" 
# topct_data = np.loadtxt("topctTestRun_topct_flow_config00000.dat",skiprows=3)
topct_data = np.loadtxt(topct_file_path,skiprows=3)

topc_file_path = "../../output/topctTestRun/flow_observables/topc/topctTestRun_topc_flow_config00000.dat" 
# topc_data = np.loadtxt("topctTestRun_topc_flow_config00000.dat",skiprows=3)
topc_data = np.loadtxt(topc_file_path,skiprows=3)

t = topct_data[:,0]
topct = topct_data[:,1:]

t_topc = topc_data[:,0]
topc = topc_data[:,1]

topct_summed = np.sum(topct,axis=1)

if np.max(np.abs(topct_summed-topc)) < 1e-16:
	print "Success!"
else:
	print "Uh-oh! Max difference is: ", np.max(np.abs(topct_summed-topc))

plt.plot(t_topc,topc,'x-',alpha=0.5,label=r"topc")
plt.plot(t,topct_summed,'-o',alpha=0.5,label=r"topct")
plt.xlabel(r"$t$")
plt.ylabel(r"$Q$")
plt.grid()
plt.legend()
plt.show()