import numpy as np, os, sys, argparse

from folderreadingtools import GetDirectoryTree, GetFolderContents

def main(batch_name,output_folder,create_setting_file=False,dryrun=False,verbose=False):
	if dryrun:
		print "%%%% dryrun %%%%"

	directory_tree = GetDirectoryTree(batch_name,output_folder,dryrun=dryrun)
	for obs in directory_tree.getFoundObservables():
		files = GetFolderContents(directory_tree,obs,flow=True)
		files.create_settings_file(dryrun=dryrun,verbose=verbose)
		files.create_perflow_data(dryrun=dryrun,verbose=verbose)

if __name__ == '__main__':
	description_string = """Small program for creating data files in per flow time format, as opposed to per configuration format."""

	parser = argparse.ArgumentParser(prog='Per flow time creator', description=description_string)

	######## Program basics #########
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	
	######## Program options ########
	parser.add_argument('batch_name', 					type=str, help='Batch name.')
	parser.add_argument('output_folder', 				type=str, help='Output folder name.')
	parser.add_argument('-csf','--create_setting_file',	default=False, action='store_true', help='Will create a seperate settings file storing info of the run.')
	parser.add_argument('--dryrun',						default=False, action='store_true', help='Dryrun where no files will be created or moved.')
	parser.add_argument('-v','--verbose',				default=False, action='store_true', help='More verbose output.')
	
	if len(sys.argv) == 1:
		args = parser.parse_args(["prodRunBeta6_1","../output/","-csf"])
	else:
		args = parser.parse_args()

	main(args.batch_name, args.output_folder, create_setting_file = args.create_setting_file, dryrun = args.dryrun, verbose = args.verbose)