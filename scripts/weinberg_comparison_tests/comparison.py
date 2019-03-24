#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
import os

def get_chroma_data(data, N=32, NT=64, NF=1001):
	"""
	Args:
		data: numpy array, chroma data
	Returns:
		w: numpy array, weinberg
		wt: numpy array, weinberg in euclidean time
	"""
	w = np.zeros(NF)
	wt = np.zeros((NF, NT))

	for i in xrange(0, NF):
		w[i] = data[i*(NT+1)]
		wt[i] = data[i*(NT+1)+1: i*(NT+1)+NT+1]

	difference = np.abs(np.sum(wt, axis=1) - w)
	print "Maximum difference between chroma summed wt and w: ", \
		np.max(difference)
	assert np.all(difference < 1e-13), "Chroma summed wt do not match w"

	return w, wt

def compare(f1, f2, N=32, NT=64, NF=1001):
	"""
	Args:
		f1: str, filename of Chroma data
		f2: str, filename of GluonAction
	"""
	f1_data = np.loadtxt(f1, usecols=[1])

	w_chroma, wt_chroma = get_chroma_data(f1_data, N=N, NT=NT, NF=NF)

	ftime, w_gluon = f2_data = np.loadtxt(f2, skiprows=3, unpack=True)

	return ftime, w_gluon, w_chroma, wt_chroma

def plot_runs(ftime, w_gluon, w_chroma, show_plot=False):
	"""
	Plots data from runs.
	"""
	plt.plot(ftime, w_gluon, label="GluonAction")
	plt.plot(ftime, w_chroma, label="Chroma")
	plt.legend()
	plt.title("Weinberg data from Chroma and GluonAction")
	plt.xlabel(r"Flow time $t$")
	plt.ylabel(r"Weinberg $W$")
	plt.grid(True)
	plt.savefig("weinberg_")
	if show_plot:
		plt.show()

def main(args):
	if len(args) == 0:
		run_num = 18
		chroma_file = "cfg3050_Qt_Ord2_Test.txt"
		gluonAction_file = "weinberg_test_data/weinbergTest%d_weinberg_flow_config00000.dat" % run_num
	else:
		assert len(args)==2, "please provide two Weinberg data files to compare."
		chroma_file = str(args[0])
		gluonAction_file = str(args[1])

	assert os.path.isfile(chroma_file), "file %s does not exist" % chroma_file
	assert os.path.isfile(gluonAction_file), "file %s does not exist" % gluonAction_file

	ftime, w_gluon, w_chroma, wt_chroma = compare(chroma_file, gluonAction_file)
	w_ratio = w_gluon / w_chroma
	w_diff = w_gluon - w_chroma

	break_point = 20

	print "i   wg             wc             wratio       wdiff"
	for i, vals in enumerate(zip(w_gluon, w_chroma, w_ratio, w_diff)):
		wg, wc, wratio, w_diff = vals
		print "%-3d %-12.8f  %-13.8f  %-12.8f %-14.8f" % (i, wg, wc, wratio, w_diff)
		if i == break_point:
			break

	plot_runs(ftime, w_gluon, w_chroma)

if __name__ == '__main__':
	import sys
	main(sys.argv[1:])