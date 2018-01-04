import os, sys, argparse, shutil, re

def natural_sort(l):
    # Natural sorting
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)',key)]
    return sorted(l,key=alphanum_key)

def main(folders,starting_integer,dryrun=False,new_batch_folder):
	# Checks all folders given are actual folders.
	CURRENT_PATH = os.getcwd()
	
	for folder in folders:
		if not os.path.isdir(folder):
			raise IOError("%s is not a valid folder." % folder)

	total_files = 0
	RENAME_FOLDERS = folders
	
	print CURRENT_PATH

	for FOLDER in RENAME_FOLDERS:

		RENAME_FOLDER_PATH = os.path.join(CURRENT_PATH,FOLDER)
		RENAME_FOLDER_FILES = natural_sort(os.listdir(RENAME_FOLDER_PATH))

		for i,file_name in enumerate(RENAME_FOLDER_FILES):
			file_path = os.path.join(CURRENT_PATH,FOLDER,file_name)
			file_base, temp = file_name.split("config")
			cfg_number, extension = temp.split(".")
			new_file_name = file_base + "config" + str(i+starting_integer) + "." + extension
			new_file_path = os.path.join(CURRENT_PATH,FOLDER,new_file_name)
			if dryrun:
				print "> mv %s %s" % (file_name,new_file_name)
			else:
				shutil.move(file_path,new_file_path)

		total_files += len(RENAME_FOLDER_FILES)
		print "Renaming %d files complete in folder %s complete." % (len(RENAME_FOLDER_FILES),FOLDER)

	print "Completed renaming a total of %d files in %d folder%s." % (total_files,len(RENAME_FOLDERS),"" if len(RENAME_FOLDERS)==1 else "s")

if __name__ == '__main__':
	# main(sys.argv)
	description_string = """
Small program intended for solving name conflicts and renaming observables if they have been created on seperate instances.
Assumes name scheme of [batch_name]beta[beta_value]_[observable_type]_[optional:flow]_config[config_number].bin
"""
	parser = argparse.ArgumentParser(prog='Observables renamer', description=description_string)

	######## Program basics #########
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	parser.add_argument('--dryrun', default=False, action='store_true', help='Dryrun to no perform any critical actions.')
	
	######## Program options ########
	parser.add_argument('folders', 					type=str, nargs='+', help='Folders to rename')
	parser.add_argument('--new_batch_folder_name',	type=str, default=False,help='Renames the batch folder to new, provided folder name.')
	parser.add_argument('starting_integer', 		type=int, help='Config integer to start new config counting from')

	args = parser.parse_args()
	# args = parser.parse_args(["--dryrun","../output/prodRunBeta6_1_15obs/flow_observables/energy","../output/prodRunBeta6_1_15obs/flow_observables/plaq","../output/prodRunBeta6_1_15obs/flow_observables/topc","485"])

	# Retrieves dryrun bool
	dryrun = args.dryrun

	main(args.folders,args.starting_integer,dryrun)
