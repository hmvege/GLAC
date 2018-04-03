import sys, os, re, struct, time

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
	bytes_to_convert = 8

	pre_time = time.clock()

	for i,fname1 in enumerate(files):
		if os.path.splitext(fname1)[-1] == ".dat": # For cases with data files
			with open(os.path.join(folder_name,fname1),"r") as file1:
				for j,fname2 in enumerate(files[i+1:]):
					# Sanity check
					if fname1 == fname2:
						sys.exit("File name duplicate: %s and %s" % (fname1, fname2))

					with open(os.path.join(folder_name,fname2),"r") as file2:
						# print "Checking %s with %s" % (fname1,fname2)
						for index, line in enumerate(zip(file1,file2)):
								if index < 3:
									continue
								if float(line[0].split(" ")[-1]) == float(line[1].split(" ")[-1]):
									equal_files.append([file1,file2])
									break
		
		elif os.path.splitext(fname1)[-1] == ".bin": # For cases with binary files
			with open(os.path.join(folder_name,fname1),"rb") as file1:
				byte1 = file1.read(bytes_to_convert)
				for j,fname2 in enumerate(files[i+1:]):
					with open(os.path.join(folder_name,fname2),"rb") as file2: 
						byte2 = file2.read(bytes_to_convert)
						if struct.unpack("d",byte1)[0] == struct.unpack("d",byte2)[0]:
							equal_files.append([file1,file2])
							continue

	post_time = time.clock()
	if len(equal_files) != 0:
		print "%d equal files found" % len(equal_files)
		for pair in equal_files:
			print pair
	else:
		print "No equal files found in folder %s\nTime used: %f seconds" % (folder_name,(post_time-pre_time))

if __name__ == '__main__':
	if len(sys.argv) == 1:
		# folder_list = [	"../data/beta6_0/flow_observables/plaq",
		# 				"../data/beta6_0/flow_observables/topc",
		# 				"../data/beta6_0/flow_observables/energy",
		# 				"../data/beta6_1/flow_observables/plaq",
		# 				"../data/beta6_1/flow_observables/topc",
		# 				"../data/beta6_1/flow_observables/energy"]
		# 				# "../data/beta6_2/flow_observables/plaq",
		# 				# "../data/beta6_2/flow_observables/topc",
		# 				# "../data/beta6_2/flow_observables/energy"]

		folder_list = [ "output/prodRunBeta6_0/flow_observables/plaq",
						"output/prodRunBeta6_0/flow_observables/topc",
						"output/prodRunBeta6_0/flow_observables/energy",
						"output/prodRunBeta6_1/flow_observables/plaq",
						"output/prodRunBeta6_1/flow_observables/topc",
						"output/prodRunBeta6_1/flow_observables/energy",]
	else:
		folder_list = sys.argv[1:]

	for folder in folder_list:
		main(folder)