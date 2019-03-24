import numpy as np, matplotlib.pyplot as plt


def getData(arr):
	tflow = arr[0::3] * 0.01
	plaq_flow = arr[1::3]
	topC_flow = arr[2::3]
	return  tflow, plaq_flow, topC_flow

# ax1 = fig1.add_subplot(311)
# t_morningstar, pf_morningstar = getFlowPlaq(plaq_morningstar_mac)
# ax1.plot(t_morningstar, pf_morningstar,"-", label="Flowed Plaquette")
# ax1.grid(True)
# # ax.set_xlabel(r"Flow time $\tau$")
# ax1.set_ylabel(r"$P_{Morningstar}$")
# ax1.set_title("Flowed plaquette value with different SU3 exponentiating methods")
# ax1.tick_params(axis='x',which='major',labelsize=8)

# t,plaq,topc = getData(results)
dat = np.loadtxt("../flowObs1.dat",dtype=float)
t = dat[:,0]
plaq = dat[:,1]
topc = dat[:,2]


fig1 = plt.figure()
ax1 = fig1.add_subplot(311)
ax1.plot(t,plaq)
ax1.set_ylabel(r"$P_t$")
ax1.set_title("Flowed plaquette and topological charge")

ax2 = fig1.add_subplot(312)
ax2.plot(t,topc,label=r"$Q = a^4\sum_x q_L^{clov}(x)$")
ax2.set_ylabel(r"$Q$")
ax2.set_xlabel(r"$t$")
ax2.legend()

ax3 = fig1.add_subplot(313)
ax3.plot(t,topc**2,label=r"$Q^2$")
ax3.set_ylabel(r"$Q^2$")
ax3.set_xlabel(r"$t$")
ax3.legend()

fig1.savefig("../figures/topological_charge_test.png",dpi=300)
plt.show()