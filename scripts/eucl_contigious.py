import numpy as np

NObs = 10 # no obs
NT = 8 # total euclidean lattice size

non_cont_mat = np.zeros((NObs, NT))

for iObs in xrange(NObs):
	for it in xrange(NT):
		non_cont_mat[iObs,it] = iObs + it

cont_mat = np.zeros(NObs*NT)

def index(_iObs, _it):
	return _iObs * NT + _it

for iObs in xrange(NObs):
	for it in xrange(NT):
		cont_mat[index(iObs, it)] = iObs + it


for iObs in xrange(NObs):
	for it in xrange(NT):
		print "%3d " % non_cont_mat[iObs, it],
	print

for iObs in xrange(NObs):
	for it in xrange(NT):
		print "%3d " % cont_mat[index(iObs, it)],
	print