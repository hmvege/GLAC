import numpy as np
import matplotlib.pyplot as plt

# Chroma data
dc = np.loadtxt("cfg3050_Qt_Ord2_Test.txt", usecols=(1))

N = 32
NT = 64
NF = 1001

w = np.zeros(NF)
wt = np.zeros((NF, NT))

for i in xrange(0, NF):
	w[i] = dc[i*(NT+1)]
	wt[i] = dc[i*(NT+1)+1: i*(NT+1)+NT+1]


difference = np.abs(np.sum(wt, axis=1) - w)
print np.max(difference)
assert np.all(difference < 1e-13), "summed wt do not match w"

# Data I generated
fname_me = "weinbergTest10_weinberg_flow_config00000.dat"
# fname_me = "weinbergTest9/flow_observables/weinberg/weinbergTest9_weinberg_flow_config00000.dat"
# fname_me = "weinbergTest6_weinberg_flow_config00000.dat"
flow_time, w_me = np.loadtxt(fname_me, skiprows=3, unpack=True)



w_ratio = w / w_me

print w
print w_me
print w_ratio
# exit(1)

plt.plot(flow_time, w_me, label="me")
plt.plot(flow_time, w, label="chroma")
plt.legend()
plt.title("Weinberg ratio between me and chroma")
plt.xlabel(r"Flow time $t$")
plt.ylabel(r"Weinberg $W$")
plt.grid(True)

# plt.figure()
# plt.plot(w_ratio)
plt.show()