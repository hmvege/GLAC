import matplotlib.pyplot as plt, numpy as np

def retrieveData(filename):
	"""
	Function for retrieving data.
	"""
	data = np.loadtxt(filename,skiprows=7)

	run_info = {}
	file = open(filename)
	for line, i in zip(file, xrange(7)):
		l = line.split()
		run_info[l[0]] = float(l[1])
	file.close()

	return data, run_info

# Data retrieval
hotStartRandomFile = "output/hotStartDataRandom.dat"
hotStartRandomData, hotStartRandomInfo = retrieveData(hotStartRandomFile)

hotStartRSTFile = "output/hotStartDataRST.dat"
hotStartRSTData, hotStartRSTInfo = retrieveData(hotStartRSTFile)

coldStartFile = "output/coldStartData.dat"
coldStartData, coldStartInfo = retrieveData(coldStartFile)

# Plotting
plt.plot(hotStartRandomData, "o", label="Hot Start(true random)", markeredgecolor="r", markerfacecolor="None")
plt.plot(hotStartRSTData, "o", label="Hot Start(RST random)", markeredgecolor="g", markerfacecolor="None")
plt.plot(coldStartData, "o", label="Cold Start", markeredgecolor="b", markerfacecolor="None")
plt.ylabel("Avg. Plaquette")
plt.xlabel("MC Updates")
plt.title("Plaquette MC History")
plt.legend(loc="best",numpoints=1)
plt.savefig("figures/MCEvolutionNonParalell.png",dpi=300)
plt.show()