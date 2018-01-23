import sys, os, re

def natural_sort(l):
	# Natural sorting
	convert = lambda text: int(text) if text.isdigit() else text.lower()
	alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)',key)]
	return sorted(l,key=alphanum_key)

def main(folder_name):
	if not os.path.isdir(folder_name):
		print "No folder name output at location %s" % folder_name
		return 0

	print "Checking files in folder %s" % folder_name

	files = natural_sort([f for f in os.listdir(folder_name) if f[0] != "."])

	equal_files = []

	number_of_double_to_check = 10

	for i,fname1 in enumerate(files):
		with open(os.path.join(folder_name,fname1),"r") as file1:
			for j,fname2 in enumerate(files[i+1:]):
				# Sanity check
				if fname1 == fname2:
					sys.exit("File name duplicate: %s and %s" % (fname1, fname2))

				with open(os.path.join(folder_name,fname2),"r") as file2:
					# print "Checking %s with %s" % (fname1,fname2)
					for index, line in enumerate(zip(file1,file2)):
						if os.path.splitext(fname1)[-1] == ".dat": # For cases with data files
							if index < 3:
								continue
							if float(line[0].split(" ")[-1]) == float(line[1].split(" ")[-1]):
								equal_files.append([file1,file2])
								break
						elif os.path.splitext(fname1)[-1] == ".bin": # For cases with binary files
							if float(line[0][0]) == float(line[1][0]):
								equal_files.append([file1,file2])
								break

	if len(equal_files) != 0:
		print "%d equal files found" % len(equal_files)
		for pair in equal_files:
			print pair
	else:
		print "No equal files found in folder %s" % folder_name

if __name__ == '__main__':
	if len(sys.argv) == 1:
		folder_list = [	"../data/beta6_0/flow_observables/plaq",
						"../data/beta6_0/flow_observables/topc",
						"../data/beta6_0/flow_observables/energy",
						"../data/beta6_1/flow_observables/plaq",
						"../data/beta6_1/flow_observables/topc",
						"../data/beta6_1/flow_observables/energy"]
						# "../data/beta6_2/flow_observables/plaq",
						# "../data/beta6_2/flow_observables/topc",
						# "../data/beta6_2/flow_observables/energy"]
	else:
		folder_list = sys.argv[1:]

	for folder in folder_list:
		main(folder)