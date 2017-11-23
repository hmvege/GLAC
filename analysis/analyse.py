from LQCDAnalyser import Bootstrap as bs, Jackknife as jk, Autocorrelation as ac
import os, numpy as np, matplotlib as plt, sys, pandas as pd


class BatchAnalyser:
	def __init__(self, output_folder, input_folder):
		self.output_folder = output_folder
		self.input_folder = input_folder

	def retrieve_data(self, run_name):
		self.run_name = run_name
		# self.base_path = os.path.join(os.path.dirname(__file__),'..')
		self.base_path = os.path.join('..')
		self.output_folder_path = os.path.join(self.base_path,self.output_folder,self.run_name)
		self.input_folder_path = os.path.join(self.base_path,self.input_folder,self.run_name)
		self.flow_observables_path = os.path.join(self.output_folder_path,'flow_observables')
		self.observables_path = os.path.join(self.output_folder_path,'observables')
		self.observables_file_path_plaq = os.path.join(self.observables_path,'plaq' + '_' + self.run_name + '.dat')
		self.observables_file_path_topc = os.path.join(self.observables_path,'topc' + '_' + self.run_name + '.dat')
		self.observables_file_path_energy = os.path.join(self.observables_path,'energy' + '_' + self.run_name + '.dat')
		self.flow_observables_path_plaq = os.path.join(self.flow_observables_path,'plaq')
		self.flow_observables_path_topc = os.path.join(self.flow_observables_path,'topc')
		self.flow_observables_path_energy = os.path.join(self.flow_observables_path,'energy')

		# Checks that folders exist
		self._check_folder_path(self.output_folder_path,True)
		self._check_folder_path(self.input_folder_path,True)

		# Checks if configuration folders exist
		self.analyse_observables = self._check_folder_path(self.observables_path)
		if self.analyse_observables:
			self._check_folder_path(self.observables_path)
			if self._check_folder_path(self.observables_file_path_plaq): self._retrieve_observables("plaq")
			if self._check_folder_path(self.observables_file_path_topc): self._retrieve_observables("topc")
			if self._check_folder_path(self.observables_file_path_energy): self._retrieve_observables("energy")

		# Checks if flow configurations folders exist
		self.analyse_flow_observables = self._check_folder_path(self.flow_observables_path)
		if self.analyse_flow_observables:
			self._check_folder_path(self.flow_observables_path)
			if self._check_folder_path(self.flow_observables_path_plaq): self._retrieve_observables(True,"plaq")
			if self._check_folder_path(self.flow_observables_path_topc): self._retrieve_observables(True,"topc")
			if self._check_folder_path(self.flow_observables_path_energy): self._retrieve_observables(True,"energy")

		# Checks if no observable-folders exist, and then promptly exists if so.
		if not self.analyse_observables and not self.analyse_flow_observables:
			raise EnvironmentError("No configuration or flow configuration folders found.")

	def _check_folder_path(self, folder_path, critical=False):
		if not os.path.isdir(folder_path):
			if critical:
				raise EnvironmentError("%s folder not found" % folder_path)
			else:
				return False
		else:
			return True

	def _retrieve_data(self, obs_type, flow=False):
		if not flow:
			self.regular_meta_data_retrieved = False
			print self.observables_file_path_plaq
			print os.listdir(self.observables_file_path_plaq)
			fname = self.observables_file_path_plaq
			print fname
		else:
			self.flow_meta_data_retrieved = False
			directory_list = os.listdir(self.flow_observables_path_plaq)
			self.NCfgs = len(directory_list)
			for i,file in enumerate(directory_list):
				# Retrieves meta-data from only one file
				if not self.flow_meta_data_retrieved:
					with open(os.path.join(self.flow_observables_path_plaq,file)) as f:
						self.beta = float((f.readline()).split()[-1])
						self.NFlows = int((f.readline()).split()[-1]) + 1 # Plus one, as we always count zeroth flow time
						self.flowEpsilon = float((f.readline()).split()[-1])
					if obs_type == "plaq":
						self.flow_plaq = np.zeros((self.NCfgs,self.NFlows))
						self.flow_time = np.zeros(self.NFlows)
						self.flow_meta_data_retrieved = True
					elif obs_type == "topc":
						self.flow_topc = np.zeros((self.NCfgs,self.NFlows))
						self.flow_time = np.zeros(self.NFlows)
						self.flow_meta_data_retrieved = True

				# Loads the columns of only one of the files
				# data_temp = np.loadtxt(os.path.join(self.flow_observables_path_plaq,file))
				if obs_type == "plaquette":
					self.flow_time, self.flow_plaq[i] = np.loadtxt(os.path.join(self.flow_observables_path_plaq,file),skiprows=3,usecols=[1,2],unpack=True)
				print self.flow_plaq

	def _retrieve_topcharge(self, flow=False):
		if flow:
			fname = self.flow_observables_path_topc
		else:
			fname = self.observables_path_topc
		data = np.loadtxt(fname)
		print fname
	

	def _retrieve_energy(self, flow=False):
		if flow:
			fname = self.observables_path_energy
		else:
			fname = self.flow_observables_path_energy
		data = np.loadtxt(fname)
		print fname
		


def main(args):
	args = ['writeLoadTest']
	BA = BatchAnalyser('output','input')
	BA.retrieve_data(args[0])

if __name__ == '__main__':
	main(sys.argv[1:])
