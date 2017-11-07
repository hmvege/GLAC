import itertools as it, numpy as np

permutations = []
for i,perm in enumerate(it.permutations('0123',4)):
	permutations.append(perm)
	# print "%2d: (%1s %1s %1s %1s)" % (i,perm[0],perm[1],perm[2],perm[3])
permutations = [map(int,perm) for perm in permutations]

permutations2 = []

c1=0
c2=0
s1=1
s2=1

for i in range(4):
	for j in range(4):
		if j == i: 
			c1+=1
			continue
		for k in range(4):
			if k == i or k == j:
				c2+=1
				continue
			for l in range(4):
				if l == k:
					c2+=1
					continue
				if l == i or l == j:
					continue
				print "(%d,%d,%d,%d) | (%d,%d)" % (i,j,k,l, 4*i + j - c1, 4*k + l - c2)
				permutations2.append([i,j,k,l])
		c2=0

if permutations == permutations2: print "TRUE: Itermethod is equal with manual method."

for perm in permutations2:
	s = 1
	for i in range(0,4):
		for j in range(i+1,4):
			s*=np.sign(perm[j]-perm[i])
	print perm, s

print "Number of permutations: %d" % len(permutations2)

l1 = []

c=0
print "Old contigious matrix allocation\ni  j  Index"
for i in range(4):
	for j in range(4):
		if i == j: 
			c += 1
			continue
		print "%-2d %-2d %-2d" % (i,j,4*i+j - c)
		l1.append(4*i+j - c)

l2 = []
print "With symmetry utilization\ni  j  Index"
c=0
for i in range(4):
	for j in range(i,4):
		if i == j: 
			c += 1
			continue
		print "%-2d %-2d %-2d" % (i,j,4*i+j-c)
		l2.append(4*i+j-c)
		# print "%-2d %-2d %-2d" % (j,i,4*j+i-c)

c=0
for j in range(4):
	for i in range(0,j+1):
		if i == j:
			c += 1
			continue
		print "%-2d %-2d %-2d" % (j,i,4*j+i-c)
		l2.append(4*j+i-c)

for i1 in l1:
	if not i1 in l2:
		print "ERROR!"

