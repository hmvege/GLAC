import os, sys, argparse, shutil, re

def natural_sort(l):
    # Natural sorting
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)',key)]
    return sorted(l,key=alphanum_key)

def main(folders,starting_integer,dryrun=False,NEW_BATCH_FOLDER=False,verbose=False,OLD_CURRENT_PATH=os.getcwd(),NEW_BATCH_NAME=False):
	# Checks all folders given are actual folders.	
	for folder in folders:
		if not os.path.isdir(folder):
			raise IOError("%s is not a valid folder." % folder)

	total_files = 0
	RENAME_FOLDERS = folders

	# Checks that the new batch folder actually exists
	if NEW_BATCH_FOLDER and not os.path.isdir(NEW_BATCH_FOLDER):
		raise IOError("%s is not a valid new batch folder." % NEW_BATCH_FOLDER)

	# Normalizes paths
	RENAME_FOLDERS = [os.path.normpath(p) for p in RENAME_FOLDERS]
	OLD_CURRENT_PATH = os.path.normpath(OLD_CURRENT_PATH)
	CURRENT_PATH = os.getcwd()
	if NEW_BATCH_FOLDER: NEW_BATCH_FOLDER = os.path.normpath(NEW_BATCH_FOLDER)

	for FOLDER in RENAME_FOLDERS:
		RENAME_FOLDER_PATH = os.path.join(OLD_CURRENT_PATH,FOLDER)

		# Sorts the folder files
		RENAME_FOLDER_FILES = natural_sort(os.listdir(RENAME_FOLDER_PATH))

		# Raises error if not all files in folder observable files of type .dat 
		if not sum([True if f.split(".")[-1] == "dat" else False for f in RENAME_FOLDER_FILES]) == len(RENAME_FOLDER_FILES):
			raise IOError("Number of files in folder %d is not of observable type .dat" % len(RENAME_FOLDER_FILES))

		# Check if a new batch folder is provided, will move observable to this
		if NEW_BATCH_FOLDER != False:
			NEW_FOLDER = NEW_BATCH_FOLDER
		else:
			NEW_FOLDER = FOLDER

		# Checks that the observable folder exists
		if NEW_BATCH_FOLDER and not os.path.isdir(NEW_FOLDER):
			raise IOError('Observable folder %s does not exist.' % NEW_FOLDER)

		# Renames all observable files in folder to provided config number
		for i,file_name in enumerate(RENAME_FOLDER_FILES):
			file_path = os.path.join(OLD_CURRENT_PATH,FOLDER,file_name)
			file_base, temp = file_name.split("config")

			# If a new batch name is provided, will replace old batch name
			if NEW_BATCH_NAME != False:
				# Finds the observable we are looking at
				observables = ["plaq","topc","energy"]
				for obs in observables:
					if obs in file_base:
						# Picks out the batch name
						batch_name, end = file_base.split("_%s" % obs)

						# Renames the file with a new batch name, and stiches the file base together
						file_base = NEW_BATCH_NAME + "_" + obs + end

						# Breaks when found the matching observable
						break
				else:
					# If observable is not found, will exit program
					raise KeyError("Observable in file %s not found among the standard observables: %s" % (file_name,", ".join(observables)))

			# Replaces the config number with a new starting integer and ensures it is on 00000 format
			cfg_number, extension = temp.split(".")
			new_file_name = file_base + "config" + "{0:0>5d}".format(i+starting_integer) + "." + extension

			# Sets up the new file path
			new_file_path = os.path.join(CURRENT_PATH,NEW_FOLDER,new_file_name)

			if os.path.isfile(new_file_path):
				raise IOError("WARNING: %s already exists" % new_file_path)

			# Performs the move/renaming of observable
			if verbose:
				print "> mv %s %s" % (file_path,new_file_path)
			else:
				print "> mv %s %s" % (os.path.relpath(file_path,OLD_CURRENT_PATH),os.path.relpath(new_file_path,CURRENT_PATH))

			# Dryrun, will not perform the move
			if not dryrun:
				shutil.move(file_path,new_file_path)

		total_files += len(RENAME_FOLDER_FILES)

		# Prints completion message depending on if have provided a new batch folder.
		if not NEW_BATCH_FOLDER:
			print "Renaming %d files in folder %s complete." % (len(RENAME_FOLDER_FILES),FOLDER)
		else:
			print "Moving %d files from folder %s to %s complete." % (len(RENAME_FOLDER_FILES),FOLDER,NEW_BATCH_FOLDER)

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
	parser.add_argument('-v', '--verbose', default=False, action='store_true', help='More verbose output.')
	
	######## Program options ########
	parser.add_argument('folders', 			type=str, nargs='+', help='Folders to rename')
	parser.add_argument('starting_integer',	type=int, help='Config integer to start new config counting from')
	parser.add_argument('--move_to_folder',	type=str, default=False,help='Moves observables to the a new folder.')
	parser.add_argument('--old_base_path',	type=str, default=os.getcwd(),help='Base folder')
	parser.add_argument('--new_batch_name',	type=str, default=False, help='Renaming files to this batch')

	args = parser.parse_args()
	# args = parser.parse_args(["--dryrun","../output/prodRunBeta6_1_15obs/flow_observables/energy","../output/prodRunBeta6_1_15obs/flow_observables/plaq","../output/prodRunBeta6_1_15obs/flow_observables/topc","485","--move_to_folder","../output/prodRunBeta6_1/flow_observables","-v"])

	main(	args.folders,
			args.starting_integer,
			dryrun=args.dryrun,
			NEW_BATCH_FOLDER=args.move_to_folder,
			verbose=args.verbose,
			OLD_CURRENT_PATH=args.old_base_path,
			NEW_BATCH_NAME=args.new_batch_name)
