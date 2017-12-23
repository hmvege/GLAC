import numpy as np, matplotlib.pyplot as plt, os, sys

class FolderTree:
	def __init__(self,flow_output_folder):
		self.flow_output_folder = flow_output_folder
		self.flow_observable_folders = os.listdir(flow_output_folder)
		self.topc_output_folder = False
		self.energy_output_folder = False
		self.plaq_output_folder = False
		if "topc" in self.flow_observable_folders: 
			self.topc_output_folder = os.path.join(flow_output_folder,"topc")
			self.topc_files = os.listdir(self.topc_output_folder)
		if "energy" in self.flow_observable_folders: 
			self.energy_output_folder = os.path.join(flow_output_folder,"energy")
			self.energy_files = os.listdir(self.energy_output_folder)
		if "plaq" in self.flow_observable_folders: 
			self.plaq_output_folder = os.path.join(flow_output_folder,"plaq")
			self.plaq_files = os.listdir(self.plaq_output_folder)

class Analyse:
	def __init__(self,folder_path,files):
		data = 


def main(argv):
	FOLDERS = FolderTree(argv[0])



if __name__ == '__main__':
	if len(sys.argv) > 1:
		args = sys.argv[1:]
	else:
		args = ["../output/prodRunBeta6_1/flow_observables/"]
	main(args)