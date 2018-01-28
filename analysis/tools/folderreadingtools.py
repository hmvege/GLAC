import sys, os, numpy as np, re

class GetDirectoryTree:
	def __init__(self,batch_name,output_folder="output",dryrun=False):
		self.flow_tree = {}
		self.obs_tree = {}
		self.CURRENT_FOLDER = os.getcwd()
		self.output_folder = output_folder
		self.observables_list = ["plaq","topc","energy"]
		self.dryrun = dryrun

		# Checks that the output folder actually exist
		if not os.path.isdir(os.path.join("..",self.output_folder)):
			raise EnvironmentError("No folder name output at location %s" % os.path.join("..",self.output_folder))
		# Retrieves folders and subfolders
		self.batch_folder = os.path.join("..",self.output_folder,batch_name)

		# Gets the regular configuration observables
		self.observables_folders = False
		obs_path = os.path.join(self.batch_folder,"observables")
		if os.path.isdir(obs_path) and len(os.listdir(obs_path)) != 0:
			self.observables_folder = obs_path
			for obs,file_name in zip(self.observables_list,os.listdir(self.observables_folder)):
				obs_path = os.path.join(self.observables_folder,file_name)
				if os.path.isfile(obs_path):
					self.obs_tree[obs] = obs_path

		# Gets paths to flow observable
		# Checks that there exists a flow observable folder
		if os.path.isdir(os.path.join(self.batch_folder,"flow_observables")):
			# Creates the flow observables path
			flow_path = os.path.join(self.batch_folder,"flow_observables")
			# Goes through the flow observables
			for flow_obs in self.observables_list:
				# Creates flow observables path
				obs_path = os.path.join(flow_path,flow_obs)
				# Checks if the flow observable path exists
				if os.path.isdir(obs_path):
					# Finds and sets the observable file paths
					flow_obs_dir_list = []
					for obs_file in os.listdir(obs_path):
						flow_obs_dir_list.append(os.path.join(obs_path,obs_file))
					# Sorts list by natural sorting
					self.flow_tree[flow_obs] = self.natural_sort(flow_obs_dir_list)
					
					# #### TEMP TEMP TEMP ####
					# self.flow_tree[flow_obs] = self.flow_tree[flow_obs][50:-300]
		print "Directory tree for batch %s built." % batch_name

		# Creates figures folder
		self.figures_path = os.path.join("..","figures",batch_name)
		if not os.path.isdir(self.figures_path):
			if self.dryrun:
				print '> mkdir %s' % self.figures_path
			else:
				os.mkdir(self.figures_path)

	@staticmethod
	def natural_sort(l):
	    # Natural sorting
	    convert = lambda text: int(text) if text.isdigit() else text.lower()
	    alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)',key)]
	    return sorted(l,key=alphanum_key)

	def getFlow(self,obs):
		"""
		Retrieves flow observable files.
		"""
		if obs in self.flow_tree.keys():
			return self.flow_tree[obs]
		else:
			raise Warning("Flow observable \"%s\" was not found in possible observables: %s" % (obs,", ".join(self.flow_tree.keys())))

	def getObs(self,obs):
		"""
		Retrieves observable files.
		"""
		if obs in self.obs_tree.keys():
			return self.obs_tree[obs]
		else:
			raise Warning("Observable \"%s\" was not found in possible observables: %s" % (obs,", ".join(self.flow_tree.keys())))

	def __str__(self):
		"""
		Prints the folder structre
		"""
		return_string = "Folder structure:"
		return_string += "\n{0:<s}".format(self.batch_folder)
		return_string += "\n{0:<s}/{1:<s}".format(self.batch_folder,"observables")
		if self.observables_folders:
			for obs,file_name in zip(self.observables_list,os.listdir(self.observables_folder)):
				return_string += "\n  {0:<s}".format(os.path.join(self.observables_folder,file_name))
		flow_path = os.path.join(self.batch_folder,"flow_observables")
		if os.path.isdir(flow_path):
			return_string += "\n  {0:<s}".format(flow_path)
			for flow_obs in (self.observables_list):
				obs_path = os.path.join(flow_path,flow_obs)
				return_string += "\n    {0:<s}".format(obs_path)
				for obs_file in os.listdir(obs_path):
					return_string += "\n      {0:<s}".format(os.path.join(obs_path,obs_file))
		return return_string

class GetFolderContents:
	"""
	Retrieves folder contents and acts as a container for data and meta-data.
	"""
	def __init__(self,folder,flow=False):
		if folder == None:
			print "    No observables found in folder: %s" % folder
		else:
			read_meta_data = True
			retrieved_flow_time = False
			retrieved_indexing = False
			N_rows_to_skip = 0
			self.meta_data = {} # Assumes meta data is the same for all files in folder
			self.data_y = []
			self.data_x = False
			if flow:
				self.data_flow_time = False	
			# Ensures we handle the data as a folder
			if type(folder) != list:
				folder = [folder]

			N_files = len(folder)

			# Goes through files in folder and reads the contents into a file
			for i,file in enumerate(folder):
				# Gets the metadata
				# print "    Reading file: %s" % file
				with open(file) as f:
					while read_meta_data:
						line = f.readline().split(" ")
						if line[0].isalpha():
							self.meta_data[str(line[0])] = float(line[-1])
							N_rows_to_skip += 1
						else:
							read_meta_data = False

				# Loads the data and places it in a list
				if flow:
					x, _x, y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
					if not retrieved_flow_time:
						self.data_flow_time = _x
						retrieved_flow_time = True
				else:
					x, y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
				
				if not retrieved_indexing:
					self.data_x = x
					retrieved_indexing = True
				self.data_y.append(y)

				# Small progressbar
				sys.stdout.write("\rData retrieved: %4.1f%% done" % (100*float(i)/float(N_files)))
				sys.stdout.flush()

			self.data_y = np.asarray(self.data_y)

			# Small progressbar
			sys.stdout.write("\rData retrieved: 100.0%% done\n")

def write_data_to_file(data,batch_name,beta,observable,folder="../output/post_analysis_data"):
	"""
	Function that write data to file.
	Args:
		data 				(numpy array)	: array containing three columns with flow times, averages and then errors in each of them
		batch_name			(str)			: file name
		observable 			(str)			: observable name
		(optional) folder 	(str)			: output folder, default is ../output/analyzed_data
	Returns:
		None
	"""
	batch_size_folder = "%s_beta%s" % (batch_name,beta)
	fname = "%s_obs%s_beta%s.txt" % (batch_name,observable,beta)
	fname_path = os.path.join(folder,batch_size_folder,fname)
	np.savetxt(fname,data,fmt="%.18f",header="t {0:<s} {0:<s}_error".format(observable))
	print "Data written to %s" % fname_path

# def join_analyzed_data_files(output)

if __name__ == '__main__':
	sys.exit("Exiting module.")