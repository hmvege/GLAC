index_list = range(4)

"""
The indices wich is printed corresponds to those missing from our list.
Will have to retrieve the dimensions of the sublattice from those
"""
for i in range(4):
	for j in range(4):
		if i!=j:
			print "i=%d j=%d  " % (i,j),
			for k in range(4):
				if k != i and k != j:
					print k,
			print ""
