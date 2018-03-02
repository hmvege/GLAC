import numpy as np
import os
import pandas as pd

def check_folder(folder_name):
	# Checks that figures folder exist, and if not will create it
	if not os.path.isdir(folder_name):
		raise IOError("Folder %s does not exist")

def main(data_folder, new_data_folder, dryrun):
	# Checks that the folders exist
	check_folder(data_folder)
	check_folder(new_data_folder)

	observable = data_folder.split("/")[-1]

	# Retrieves data files
	folder_files = sorted(os.listdir(data_folder))
	new_folder_files = sorted(os.listdir(new_data_folder))

	# Sets up dictionary for transporting the meta-data
	meta_data = {}

	for i, files in enumerate(zip(folder_files,new_folder_files)):
		file, new_file = files

		N_rows_to_skip = 0
		read_meta_data = True

		file_path = os.path.join(data_folder,file)
		new_file_path = os.path.join(new_data_folder,new_file)

		# Read meta-data from topct (they are the same for all files)
		with open(file_path) as f:
			# Reads in meta data as long as the first element on the line is a string
			while read_meta_data:
				line = f.readline().split(" ")
				if line[0].isalpha():
					meta_data[str(line[0])] = float(line[-1])
					N_rows_to_skip += 1
				else:
					# Stores number of rows(if we are on old or new data reading)
					N_rows = len(line)

					NT = N_rows - 1 # Temporal dimension
					NFlows = int(meta_data["NFlows"])

					# Exits while loop
					read_meta_data = False

		# Read old folder data
		# Uses pandas to read data (quick!)
		# Sets up header names
		header_names = list("t")
		header_names[1:] = ["t%d" % i for i in range(NT)]

		data = pd.read_csv(file_path, skiprows=N_rows_to_skip, sep=" ", names=header_names, header=None)

		data_array = np.zeros((NFlows,NT))
		summed_array = np.zeros(NFlows)

		for i,itemporal in enumerate(header_names[1:]):
			data_array[:,i] = data[itemporal]

		summed_array = np.sum(data_array,axis=1)

		if not dryrun:
			with open(new_file_path,'w') as f:
				# Write new meta-data
				f.write("beta %.2f" % meta_data["beta"])
				f.write("\nNFlows %d" % NFlows)
				f.write("\nFlowEpsilon %.2f" % meta_data["FlowEpsilon"])
				
				# Write new topc data
				for tau, obs in zip(data["t"].values,summed_array):
					f.write("\n%.16f %.16f" % (tau,obs))
		else:
			print "\n\nWriting to file:\n"
			print "beta %.2f" % meta_data["beta"]
			print "NFlows %d" % meta_data["NFlows"]
			print "FlowEpsilon %.2f" % meta_data["FlowEpsilon"]
			for tau, obs in zip(data["t"].values,summed_array):
					print "%.16f %.16f" % (tau,obs)

		print "File %s written. " % new_file_path

	print "Summing and rewriting complete."

if __name__ == '__main__':
	dryrun = False

	topct_folder = "../data5/beta62/flow_observables/topct"
	new_topc_folder = "../data5/beta62/flow_observables/topc"
	main(topct_folder, new_topc_folder, dryrun)