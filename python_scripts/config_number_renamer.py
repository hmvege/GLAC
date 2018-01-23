import sys, os, re, shutil

def folder_looper(folder, dryrun, verbose = True):
	print "Renaming files in folder: %s" % folder
	for file in os.listdir(folder):
		# Puts together file path
		file_path = os.path.join(folder,file)

		# If sub folder is a folder as well, we recursively enter it. Afterwards, we exit it
		if os.path.isdir(file_path):
			folder_looper(file_path,dryrun)
			continue
		
		# Splits file name into file head and extension
		file_head, file_ext = os.path.splitext(file)

		# Error catching in case file format is unknown or similar
		if file_ext not in ['.bin','.dat']:
			print "Unrecognized file format %s for file %s" % (file_ext,file_path)
			continue
		
		# Checks that we have a file that contains a config number in its name
		if not "config" in file_head:
			print "Unrecognized file type: no configuration found in file_head %s" % file_head
			continue

		# Retrieves the old config number
		files_stripped = [c for c in re.split('(\d+)',file_head) if c.isdigit()]

		# If we already have the correct numbering counting, we skip
		if len(files_stripped[-1]) == 5:
			continue

		# Converts string to digits
		files_stripped = [int(c) for c in files_stripped]

		# Plucks out the old configuration number
		old_config_number = files_stripped[-1]

		# Sets up new configuration number
		new_config_number = "{:0>5d}".format(files_stripped[-1])

		# Sets up new file name
		new_file_name = file_head.split('config')[0] + 'config' + new_config_number + file_ext

		# Sets up new file name with folder base
		new_file_path = os.path.join(folder,new_file_name)

		# Performs file renaming, or possible the dryrun equivalent
		if dryrun:
			print "> move %s %s" % (file_path,new_file_path)
		else:
			if verbose:
				print "> move %s %s" % (file_path,new_file_path)
			shutil.move(file_path,new_file_path)

def main(folder,dryrun = False, verbose = False):
	if dryrun: print "{0:<s} Dryrun mode on {0:<s}\n".format("="*20)

	# Error checking
	if not os.path.isdir(folder):
		print "%s is not a valid folder" % folder
		return

	# Loops through folders inside recursively, and renames accordingly
	folder_looper(folder,dryrun,verbose)

if __name__ == '__main__':
	dryrun = False
	verbose = True
	if len(sys.argv) == 1:
		# folder_list = ["../output/ubuntu_test_run"]
		folder_list = ["../output/test_run_new_counting/"]
	else:
		folder_list = sys.argv[1:]

	for folder in folder_list:
		main(folder,dryrun=dryrun,verbose=verbose)