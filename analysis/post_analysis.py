import numpy as np, matplotlib.pyplot as plt, sys, os

class DataReader:
	"""
	Small class for reading post analysis data
	"""
	def __init__(self,post_analysis_folder):
		batch_folders = self._get_folder_content(post_analysis_folder)

		data = {}
		for batch in batch_folders:
			batch_folder = os.path.join(post_analysis_folder,batch)
			batch_files = [os.path.join(batch_folder,i) for i in self._get_folder_content(batch_folder)]

			# Observable data list
			observable_data = []

			for bfile in batch_files:
				# Retrieves meta data from header
				meta_data = ""
				with open(bfile) as f:
					header_content = f.readline().split(" ")[1:]
					meta_data = [h.split("\n")[0] for h in header_content]

				# Retrieves observable from meta_data
				observable = meta_data[1] 

				observable_data.append({"observable": observable, "data": np.loadtxt(bfile)})

			# Loads data and stores it in a dictionary for all the data
			data[batch] = observable_data

	def write_batch_to_single_file(self):
		None

	@staticmethod
	def _get_folder_content(folder):
		if not os.path.isdir(folder):
			raise IOError("No folder by the name %s found." % folder)
		else:
			return os.listdir(folder)

class TopSusPostAnalysis:
	None

class EnergyPostAnalysis:
	None

def main(args):
	"""
	Args should be post-analysis folder
	"""
	# Loads data from post analysis folder
	data = DataReader(args[0])
	print "Retrieving data from folder: %s" % args[0]

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