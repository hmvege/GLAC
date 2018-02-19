import sys, os, numpy as np, re, pandas as pd, time

class GetDirectoryTree:
	def __init__(self,batch_name,batch_folder,dryrun=False):
		self.flow_tree = {}
		self.obs_tree = {}
		self.CURRENT_FOLDER = os.getcwd()
		self.data_batch_folder = batch_folder
		self.observables_list = ["plaq","topc","energy"]
		self.batch_name = batch_name
		self.dryrun = dryrun

		# Checks that the output folder actually exist
		if not os.path.isdir(os.path.join("..",self.data_batch_folder)):
			raise EnvironmentError("No folder name output at location %s" % os.path.join("..",self.data_batch_folder))
		# Retrieves folders and subfolders
		self.batch_name_folder = os.path.join("..",self.data_batch_folder,batch_name)

		# Gets the regular configuration observables
		self.observables_folders = False
		obs_path = os.path.join(self.batch_name_folder,"observables")
		if os.path.isdir(obs_path) and len(os.listdir(obs_path)) != 0:
			self.observables_folder = obs_path
			for obs,file_name in zip(self.observables_list,os.listdir(self.observables_folder)):
				obs_path = os.path.join(self.observables_folder,file_name)
				if os.path.isfile(obs_path):
					self.obs_tree[obs] = obs_path

		# Gets paths to flow observable
		# Checks that there exists a flow observable folder
		if os.path.isdir(os.path.join(self.batch_name_folder,"flow_observables")):
			# Creates the flow observables path
			flow_path = os.path.join(self.batch_name_folder,"flow_observables")

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

		# Creates figures folder
		if os.path.split(self.CURRENT_FOLDER)[-1] == "tools":
			self.figures_path = os.path.join("..","..","figures",batch_name)
		elif os.path.split(self.CURRENT_FOLDER)[-1] == "analysis":
			self.figures_path = os.path.join("..","figures",batch_name)
		else:
			raise OSError("Current folder path not recognized: %s" % self.CURRENT_FOLDER)
		
		check_folder(self.figures_path,self.dryrun,verbose=True)

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

	def getFoundObservables(self):
		"""
		Returns list over all found observables.
		"""
		return self.observables_list

	def __str__(self):
		"""
		Prints the folder structre
		"""
		return_string = "Folder structure:"
		return_string += "\n{0:<s}".format(self.batch_name_folder)
		return_string += "\n  {0:<s}/{1:<s}".format(self.data_batch_folder,"observables")
		if self.observables_folders:
			for obs,file_name in zip(self.observables_list,os.listdir(self.observables_folder)):
				return_string += "\n    {0:<s}".format(os.path.join(self.observables_folder,file_name))
		flow_path = os.path.join(self.batch_name_folder,"flow_observables")
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
	def __init__(self,file_tree,observable,flow=False):
		# Retrieves file from file tree
		files = file_tree.getFlow(observable)

		# Stores the file tree as a global constant for later use by settings and perflow creator
		self.file_tree = file_tree

		# Stores the observable name
		self.observable = observable

		# t1 = time.clock()

		if files == None:
			print "    No observables of type %s found in folder: %s" % (observable,folder)
		else:
			# Bools to ensure certain actions are only done once
			read_meta_data = True
			retrieved_flow_time = False
			retrieved_indexing = False

			# Number of rows to skip after meta-data has been read
			N_rows_to_skip = 0

			# Long-term storage variables
			self.meta_data = {}
			self.data_y = []
			self.data_x = False

			# If we are not 
			if flow:
				self.data_flow_time = False	

			# Ensures we handle the data as a folder
			if type(files) != list:
				self.files = [files]
			else:
				self.files = files

			# Number of files is the length of files in the the folder
			N_files = len(self.files)

			# Goes through files in folder and reads the contents into a file
			for i,file in enumerate(self.files):
				# print "    Reading file: %s" % file

				# Gets the metadata
				with open(file) as f:
					# Reads in meta data as long as the first element on the line is a string
					while read_meta_data:
						line = f.readline().split(" ")
						if line[0].isalpha():
							self.meta_data[str(line[0])] = float(line[-1])
							N_rows_to_skip += 1
						else:
							read_meta_data = False

				# Loads the data and places it in a list
				if flow:
					# Uses numpy to load flow data (slow!)
					# x, _x, y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
					
					# Uses pandas to read data (quick!)
					data = pd.read_csv(file,skiprows=N_rows_to_skip,sep=" ",names=["t","sqrt8t",self.observable],header=None)

					# Only retrieves flow once
					if not retrieved_flow_time:
						# self.data_flow_time = _x # This is the a*sqrt(8*t), kinda useless
						self.data_flow_time = data["sqrt8t"].values # Pandas
						retrieved_flow_time = True
				else:
					# Uses numpy to load non-flown data (slow!)
					# x, y = np.loadtxt(file,skiprows=N_rows_to_skip,unpack=True)
					
					# Uses pandas to load non-flown data (quick!)
					data = pd.read_csv(file,skiprows=N_rows_to_skip,sep=" ",names=["t","sqrt8t",self.observable],header=None)
				
				# Only retrieves indexes/flow-time*1000 once
				if not retrieved_indexing:
					# self.data_x = x # Numpy
					self.data_x = data["t"].values # Pandas
					retrieved_indexing = True
				# self.data_y.append(y)
				self.data_y.append(data[self.observable].values)

				# # Small progressbar
				# sys.stdout.write("\rData retrieved: %4.1f%% done" % (100*float(i)/float(N_files)))
				# sys.stdout.flush()

			# Converts data to a numpy array
			self.data_y = np.asarray(self.data_y)

			# # Small progressbar
			# sys.stdout.write("\rData retrieved: 100.0%% done\n")

		# time_used = time.clock()-t1	
		# print "Data reading: time used with function: %.10f secs/ %.10f minutes" % (time_used, time_used/60.)	

	def create_perflow_data(self,dryrun=False,verbose=False):
		# Creating per flow folder
		per_flow_folder = os.path.join("..",self.file_tree.data_batch_folder,"perflow")
		check_folder(per_flow_folder,dryrun,verbose=verbose)

		# Creates observable per flow folder
		per_flow_observable_folder = os.path.join(per_flow_folder,self.observable)
		check_folder(per_flow_observable_folder,dryrun,verbose=verbose)

		# Retrieving number of configs and number of flows
		NConfigs,NFlows = self.data_y.shape

		# Re-storing files in a per flow format
		for iFlow in xrange(NFlows):
			# Setting up new per-flow file
			flow_file = os.path.join(per_flow_folder,self.observable,self.file_tree.batch_name + "_flow%05d.dat" % iFlow)

			# Saving re-organized data to file
			if not dryrun:
				np.savetxt(flow_file,self.data_y[:,iFlow],fmt="%.16f",header="t %f \n%s" % (iFlow*self.meta_data["FlowEpsilon"],self.observable))

			# Prints message regardless of dryrun and iff
			if verbose:
				print "    %s created." % flow_file

		print "Per flow data for observable %s created." % self.observable

	def create_settings_file(self,dryrun=False,verbose=False):
		"""
		Function for storing run info.
		"""		
		# Checking that settings file does not already exist
		setting_file_path = os.path.join(self.file_tree.batch_folder,"run_settings_%s.txt" % self.observable)

		# Creating string to be passed to info file
		info_string = ""
		info_string += "Batch %s" % self.file_tree.batch_name
		info_string += "\nObservables %s" % " ".join(self.file_tree.getFoundObservables())
		info_string += "\nNConfigs %d" % self.data_y.shape[0]
		for key in self.meta_data:
			info_string += "\n%s %s" % (key, self.meta_data[key])

		# Creating file
		if not dryrun:
			with open(setting_file_path,"w") as f:
				# Appending batch name
				f.write(info_string)

		# Prints file if we have not dryrun or if verbose option is turned on
		if verbose:
			print "\nSetting file:", info_string, "\n"
		
		print "Setting file %s created." % setting_file_path

def check_folder(folder_name, dryrun, verbose = False):
	# Checks that figures folder exist, and if not will create it
	if not os.path.isdir(folder_name):
		if dryrun or verbose:
			print "> mkdir %s" % folder_name
		if not dryrun:
			os.mkdir(folder_name)


def write_data_to_file(analysis_object, post_analysis_folder = "../output/post_analysis_data", dryrun = False):
	"""
	Function that write data to file.
	Args:
		analysis_object		(FlowAnalyser)	: object of analyser class
		(optional) folder 	(str)			: output folder, default is ../output/analyzed_data
	Returns:
		None
	"""
	# Retrieves beta value and makes it into a string
	beta_string = str(analysis_object.beta).replace(".","_")
	
	# Ensures that the post analysis data folder exists
	check_folder(post_analysis_folder,dryrun,verbose=True)

	# Sets up batch folder
	batch_folder_path = os.path.join(post_analysis_folder,analysis_object.batch_data_folder)
	check_folder(batch_folder_path,dryrun,verbose=True)

	# Sets up beta value folder
	beta_folder_path = os.path.join(batch_folder_path,"beta" + beta_string)
	check_folder(beta_folder_path,dryrun,verbose=True)

	# Retrieves analyzed data
	x = analysis_object.x*analysis_object.data.meta_data["FlowEpsilon"]
	y_org = analysis_object.unanalyzed_y
	y_err_org = analysis_object.unanalyzed_y_std*analysis_object.autocorrelation_error_correction
	y_bs = analysis_object.bs_y
	y_err_bs = analysis_object.bs_y_std*analysis_object.autocorrelation_error_correction
	y_jk = analysis_object.jk_y
	y_err_jk = analysis_object.jk_y_std*analysis_object.autocorrelation_error_correction

	# Stacks data to be written to file together
	data = 	np.stack((x, y_org, y_err_org, y_bs, y_err_bs, y_jk, y_err_jk),axis=1)
	
	# Retrieves compact analysis name 
	observable = analysis_object.observable_name_compact

	# Sets up file name and file path
	fname = "%s.txt" % observable
	fname_path = os.path.join(beta_folder_path,fname)

	# Saves data to file
	if not dryrun:
		np.savetxt(fname_path,data,fmt="%.16f",header="observable %s beta %s\nt orginal original_error bs bs_error jk jk_error" % (observable,beta_string))
	print "Data written to %s" % fname_path

def write_raw_analysis_to_file(raw_data, analysis_type, observable, data_batch_folder, beta, post_analysis_folder = "../output/post_analysis_data", dryrun = False):
	"""
	Function that writes raw analysis data to file, either bootstrapped or jackknifed data.
	Args:
		raw_data (numpy float array, NFlows x NBoot)
		analysis_type (str)
		observable (str)
		(optional) post_analysis_folder (str)
		(optional) dryrun (bool)
	"""
	# Ensures that the post analysis data folder exists
	check_folder(post_analysis_folder,dryrun,verbose=True)

	# Sets up data batch folder
	data_batch_folder_path = os.path.join(post_analysis_folder,data_batch_folder)
	check_folder(data_batch_folder_path,dryrun,verbose=True)

	# Sets up beta folder
	beta_folder_path = os.path.join(data_batch_folder_path,"beta" + str(beta).replace(".","_"))
	check_folder(beta_folder_path,dryrun,verbose=True)

	# Sets type type of analysis folder
	analysis_folder_path = os.path.join(beta_folder_path,analysis_type)
	check_folder(analysis_folder_path,dryrun,verbose=True)

	# Creates file name
	file_name = observable
	file_name_path = os.path.join(analysis_folder_path,file_name)

	# Stores data as binary output
	np.save(file_name_path,raw_data)
	print "Analysis %s for observable %s for beta %.2f stored as binary data at %s.npy" % (analysis_type,observable,beta,file_name_path)

if __name__ == '__main__':
	sys.exit("Exiting module.")