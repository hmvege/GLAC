numprocs = 523

maxRank = numprocs

binaryCounter = 1
counter = 0

while maxRank - binaryCounter != 0:
	binaryCounter = 2**counter
	counter += 1

	rest = maxRank % binaryCounter
	if rest != 0:
		maxRank -= rest
		print "Subtracting ", rest

	print binaryCounter, maxRank

