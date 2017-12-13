
for direction in range(8):
	print "\nDirection: %d \nPositive: %d \nReceiving from: %d" % (direction, direction % 2, ((direction + 1) % 2) + direction / 2 * 2)
	print "mu=%d %s" % (direction / 2, "FORWARDS" if direction % 2 == 0 else "BACKWARDS")
	print "Cube: %d %d %d" % ((direction / 2 + 1) % 4, (direction / 2 + 2) % 4, (direction / 2 + 3) % 4)

# (dir / 2 + 1) % 4