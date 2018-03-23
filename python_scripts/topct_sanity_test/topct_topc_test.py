import numpy as np, os, matplotlib.pyplot as plt

print os.getcwd(), os.listdir(os.getcwd())

# topct_file_path = "../../output/topctTestRun/flow_observables/topct/topctTestRun_topct_flow_config00000.dat" 
# topct_data = np.loadtxt("../../data5/beta61/flow_observables/topct/prodRunBeta6_1_topct_flow_config00501.dat",skiprows=3)
topct_data = np.loadtxt("../../data6/beta60/flow_observables/topct/prodRunBeta6_0_topct_flow_config00001.dat",skiprows=3)
# topct_data = np.loadtxt("../../output/lattice_field_topct_test/flow_observables/topct/lattice_field_topct_test_topct_flow_config00001.dat",skiprows=3)
# topct_data = np.loadtxt("prodRunBeta6_1_topct_flow_config00400.dat",skiprows=3)
# topct_data = np.loadtxt(topct_file_path,skiprows=3)

# topc_file_path = "../../output/topctTestRun/flow_observables/topc/topctTestRun_topc_flow_config00000.dat" 
# topc_data = np.loadtxt("../../data5/beta61/flow_observables/topc/prodRunBeta6_1_topc_flow_config00501.dat",skiprows=3)
topc_data = np.loadtxt("../../data6/beta60/flow_observables/topc/prodRunBeta6_0_topc_flow_config00001.dat",skiprows=3)
# topc_data = np.loadtxt("../../output/lattice_field_topct_test/flow_observables/topc/lattice_field_topct_test_topc_flow_config00001.dat",skiprows=3)
# topc_data = np.loadtxt("prodRunBeta6_1_topc_flow_config00400.dat",skiprows=3)
# topc_data = np.loadtxt(topc_file_path,skiprows=3)

if topc_data.shape[0] > topct_data.shape[0]:
	topc_data = topc_data[:topct_data.shape[0]]

t = topct_data[:,0]
topct = topct_data[:,1:]

t_topc = topc_data[:,0]
topc = topc_data[:,1]

t = np.sqrt(8*t)
t_topc = np.sqrt(8*t_topc)

topct_summed = np.sum(topct, axis=1)


if np.max(np.abs(topct_summed-topc)) < 1e-16:
	print "Success!"
else:
	print "Uh-oh! Max difference is: ", np.max(np.abs(topct_summed-topc))

fig = plt.figure(1)
ax1 = fig.add_subplot(211)
ax1.plot(t_topc,topc,'x-',alpha=0.5,label=r"topc")
ax1.plot(t,topct_summed,'-o',alpha=0.5,label=r"topct")
ax1.set_xlabel(r"$t$")
ax1.set_ylabel(r"$Q$")
ax1.grid()
ax1.legend()

ax2 = fig.add_subplot(212)
ax2.plot(t_topc, np.abs(topc-topct_summed), label=r"|topc-topct|")
ax2.grid()
ax2.legend()
# ax2.set_ylim(10**(-16),0)
plt.show()