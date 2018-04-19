import numpy as np

# Number of colors / N in SU(N)
Nc = 3

def index_mapper(i, j):
	return int(i!=0)*Nc + j-1

next_index = lambda _i: _i % Nc + 1

mu = 0
for nu in range(1,4):
	rho = next_index(nu)
	sigma = next_index(rho)

	print mu, nu, rho, sigma, "map1: %d" % index_mapper(mu, nu), "map2: %d" % index_mapper(rho, sigma)


