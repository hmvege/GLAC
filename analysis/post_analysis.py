import numpy as np, matplotlib.pyplot as plt, sys, os

class DataReader:
	"""
	Small class for reading post analysis data
	"""
	def __init__(self,folder):
		if not os.path.isdir(folder):
			raise IOError("No folder by the name %s found." % folder)
		else:
			data_folder = os.listdir(folder) # Recursively scrape folders

class TopSusPostAnalysis:
	None

class EnergyPostAnalysis:
	None

def main(args):
	print "Retrieving data from folder: %s" % args[0]
	data = DataReader(args[0])
	# Loads data from all sources for each of the type of observables, 

	# Rewrites all of the data to a single file for sharing with giovanni

	# Plots topsus

	# Retrofits the topsus for continiuum limit

	# Plots energy

	# Retrofits the energy for continiuum limit

if __name__ == '__main__':
	if len(sys.argv[1:]) == 1:
		args = sys.argv[1:]
	else:
		args = ["../output/post_analysis_data"]
	main(args)