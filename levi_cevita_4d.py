import itertools as it, numpy as np

permutations = []
for i,perm in enumerate(it.permutations('0123',4)):
	permutations.append(perm)
	print "%2d: (%1s %1s %1s %1s)" % (i,perm[0],perm[1],perm[2],perm[3])
permutations = [map(int,perm) for perm in permutations]

permutations2 = []
for i in range(4):
	for j in range(4):
		if j == i: continue
		for k in range(4):
			if k == i or k == j: continue
			for l in range(4):
				if l == i or l == j or l == k: continue
				permutations2.append([i,j,k,l])

if permutations == permutations2: print "TRUE: Itermethod is equal with manual method."

for perm in permutations:
	s = 1
	for i in range(0,4):
		for j in range(i+1,4):
			s*=np.sign(perm[j]-perm[i])
	print perm, s