import numpy as np
import matplotlib.pyplot as plt
dat = np.load("data8/beta62/post_analysis_data/qtq0eff/tflow0999/bootstrap/qtq0eff.npy")
print "dat.shape", dat.shape
print "dat[:,0]", dat[:,0]
print "np.roll(dat,-1,axis=0)[:,0]", np.roll(dat,-1,axis=0)[:,0]
print "np.mean(dat / np.roll(dat,-1,axis=0),axis=1):", np.mean(dat / np.roll(dat,-1,axis=0),axis=1)
# print "np.log(np.mean(dat / np.roll(dat,-1,axis=0),axis=1))", np.log(np.mean(dat / np.roll(dat,-1,axis=0),axis=1))
# print "np.log(np.mean(np.roll(dat,-1,axis=0)", np.log(np.mean(np.roll(dat,-1,axis=0)/dat,axis=1))


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
def plot_data(data_type):
	# Data gathering
	data = np.load("data8/beta62/post_analysis_data/qtq0eff/tflow0999/%s/qtq0eff.npy" % data_type)
	log_data = np.log(data/np.roll(data, -1, axis=0))
	y = np.mean(log_data, axis=1)
	yerr = np.std(log_data, axis=1)

	eucl_t = 32
	ax1.hist(data[eucl_t],label=data_type)
	ax1.legend()
	ax1.set_title(r"$t_e=%d$" % eucl_t)
	ax1.set_xlabel(r"$aM_{eff}$")
	fig1.savefig("eff_mass_hist_b62.png", dpi=400)

	# Plot commands
	ax2.errorbar(range(data.shape[0]), y, yerr=yerr, capsize=5, fmt="_", ls=":",label=data_type)
	ax2.set_title("Effective mass")
	ax2.set_ylabel(r"$aM_{eff}$")
	ax2.set_xlabel(r"$t_e/a$")
	ax2.legend()
	fig2.savefig("eff_mass_b62.png", dpi=400)

data_types = ["unanalyzed", "bootstrap", "jackknife"]

for itype in data_types:
	plot_data(itype)

plt.show()