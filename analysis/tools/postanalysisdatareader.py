from folderreadingtools import check_folder
import os, numpy as np, copy

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a # fermi

class PostAnalysisDataReader:
	"""
	Small class for reading post analysis data
	"""
	def __init__(self,post_analysis_folder,verbose=False):
		self.post_analysis_folder = post_analysis_folder
		self.verbose = verbose

		# Dictionary variable to hold all the data sorted by batches
		self.data_batches = {}

		# Dictionaries to hold the raw bootstrapped/jackknifed data and autocorrelation data
		self.bs_data_raw = {}
		self.jk_data_raw = {}
		self.ac_corrections = {}

		# Different types of analysis
		self.analysis_types = ["unanalyzed","jackknife","bootstrap"]

		# Binary folder types available
		self.binary_folder_types = ["jackknife","bootstrap","autocorrelation"]

		# Variable to store if we have retrieved flow time or not
		self.retrieved_flow_time = False

		# Number of betas variable
		self.N_betas = 0

		# Data batch name
		self.data_batch_name = os.path.split(self.post_analysis_folder)[-1]

		# Iterates over the different beta value folders
		for beta_folder in self._get_folder_content(self.post_analysis_folder):
			# Construct beta folder path
			beta_folder_path = os.path.join(self.post_analysis_folder,beta_folder)

			# Retrieves beta from folder name
			beta = float(beta_folder.strip("beta").replace("_","."))

			# Dictionary to store observable data in
			observable_data = {}

			# Loops over files and folders inside the beta folder
			for folder in self._get_folder_content(beta_folder_path):
				# Construct folder path
				folder_path = os.path.join(beta_folder_path,folder)

				# Retrieves data deepending on type
				if folder == "bootstrap":
					# Gets binary bs data
					self.bs_data_raw[beta] = self._get_bin_dict(folder_path)
				elif folder == "jackknife":
					# Gets binary jk data
					self.jk_data_raw[beta] = self._get_bin_dict(folder_path)
				elif folder == "autocorrelation":
					# Gets binary ac data
					self.ac_corrections[beta] = self._get_bin_dict(folder_path)
				else:
					# Gets observable name
					observable_name = os.path.splitext(folder)[0]
					
					# Retrieves the observable data
					observable_data[observable_name] = self._get_beta_observable_dict(folder_path,beta)

			# Stores batch data
			self.data_batches[beta] = copy.deepcopy(observable_data)

			# Add another beta value
			self.N_betas += 1

			# Frees memory
			del observable_data
		
		# Reorganizes data to more ease-of-use type of data set
		self._reorganize_data()

	def _get_beta_observable_dict(self,observable_file,beta):
		# Retrieves metadata
		# Make it so one can retrieve the key as meta_data[i] and then value as meta_data[i+1]
		meta_data = self._get_meta_data(observable_file)

		# Temporary methods for getting observable name and beta value, as this will be put into the meta data
		obs = os.path.split(os.path.splitext(observable_file)[0])[-1]

		# Dictionary to store all observable data in
		obs_data = {}

		# Loads data into temporary holder
		retrieved_data = np.loadtxt(observable_file)

		# Puts data into temporary holding facilities
		t 			= retrieved_data[:,0]
		y 			= retrieved_data[:,1]
		y_error 	= retrieved_data[:,2]
		bs_y 		= retrieved_data[:,3]
		bs_y_error 	= retrieved_data[:,4]
		jk_y 		= retrieved_data[:,5]
		jk_y_error 	= retrieved_data[:,6]


		# Stores data into dictionaries
		unanalyzed_data = {"y": y, "y_error": y_error}
		bs_data = {"y": bs_y, "y_error": bs_y_error}
		jk_data = {"y": jk_y, "y_error": jk_y_error}

		# Stores observable data
		obs_data["beta"] 		= copy.deepcopy(beta)
		obs_data["unanalyzed"] 	= copy.deepcopy(unanalyzed_data)
		obs_data["bootstrap"] 	= copy.deepcopy(bs_data)
		obs_data["jackknife"] 	= copy.deepcopy(jk_data)

		# Stores flow time in a seperate variable
		if not self.retrieved_flow_time:
			self.flow_time = copy.deepcopy(t)
			self.retrieved_flow_time = True

		if self.verbose:
			print "Data retrieved from %s" % observable_file

		# Frees memory
		del retrieved_data

		return obs_data

	def _get_bin_dict(self,folder):
		"""
		Gets binary data files
		"""
		observable_dict = {}

		# Retrieves the different bootstrapped observables
		for observable_file in self._get_folder_content(folder):
			# Gets observable name
			observable_name = os.path.splitext(observable_file)[0]

			# Populates dictionary
			observable_dict[observable_name] = np.load(os.path.join(folder,observable_file))

		return observable_dict

	@staticmethod
	def _get_meta_data(file):
		# Retrieves meta data from header or file
		meta_data = ""
		with open(file) as f:
			header_content = f.readline().split(" ")[1:]
			meta_data = [h.split("\n")[0] for h in header_content]
		return meta_data

	def _reorganize_data(self):
		# Reorganizes the data into beta-values and observables sorting
		self.data_observables = {}

		# Sets up new dictionaries by looping over batch names
		for beta in self.data_batches:
			# Loops over observable names
			for observable_name in self.data_batches[beta]:
				# Creates new sub-dictionary ordered by the observable name
				self.data_observables[observable_name] = {}

		# Places data into dictionaries
		for beta in self.data_batches:
			# Loops over the batch observable
			for observable_name in self.data_batches[beta]:
				# Stores the batch data in a sub-dictionary
				self.data_observables[observable_name][beta] = self.data_batches[beta][observable_name]


	def write_batch_to_single_file(self):
		"""
		Writes all unanalyzed data to a single file
		"""
		# Creates universal output folder
		file_folder = os.path.join(self.post_analysis_folder,"..","..","universal_output") # Should be in output
		print "="*100,"\nWriting data to a universal output format in %s" % file_folder
		check_folder(file_folder,False,True)

		# Variable for checking if we have retrieved the flow time or not
		retrieved_flow_time = False

		# Loops over the different possible analyses
		for analysis_type in self.analysis_types:
			# Loops over the batches
			for beta in self.data_batches:
				# Temporary data output for gathering data from arrays into a simple list
				data_output = []

				# Header for the resulting array
				header_output = ["t","sqrt8t"]

				# Retrieves t
				data_output.append(self.flow_time)

				# Retrieves sqrt(8t)
				data_output.append(np.sqrt(8*self.flow_time))

				for observable_name in self.data_batches[beta]:
					# Retrieves the unanalyzed observable data
					data_output.append(self.data_batches[beta][observable_name][analysis_type]["y"])
					data_output.append(self.data_batches[beta][observable_name][analysis_type]["y_error"])

					# Adds observable name and error to header
					header_output.append(observable_name)
					header_output.append(observable_name+"_err")

				# Writing array to file
				file_path = os.path.join(file_folder,str(beta).replace(".","_") + "_" + analysis_type) + ".txt"
				np.savetxt(file_path,np.asarray(data_output),fmt="%.18f",header=" ".join(header_output))
				print "Batch '%s' data written to file %s." % (beta,os.path.splitext(os.path.split(file_path)[-1])[0])

		# A single delimiter marking the end of universal output writing
		print "="*100

	@staticmethod
	def _get_folder_content(folder):
		if not os.path.isdir(folder):
			raise IOError("No folder by the name %s found." % folder)
		else:
			return [ f for f in os.listdir(folder) if not f.startswith(".") ]

if __name__ == '__main__':
	import sys
	sys.exit("Error: PostAnalysisDataReader is not a standalone program.")