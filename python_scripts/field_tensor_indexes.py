import numpy as np

# Number of colors / N in SU(N)
Nc = 3
max_int = 6

# _mapper = lambda i, j: int(i!=0)*Nc + j - 1 + (j-i-1)%Nc
_mapper = lambda i, j: int(i!=0)*Nc + j - 1
# _mapper = lambda i, j: ((j-i)%Nc + 1)/Nc  

s = {
	0 :{
		1: 0,
		2: 1,
		3: 2,
	},
	1 : {
		2: 3,
		3: 4,
	},
	2 : {
		3: 5,
	},
}

print s[2][3]
# exit()

def index_mapper(i, j):
	return _mapper(i,j)

# def index_mapper(i, j):
# 	if i < j:
# 		# return _mapper(i,j)
# 		return s[i][j]
# 	else:
# 		return s[j][i]

next_index = lambda _i: _i % Nc + 1

# for i in range(4):
# 	for j in range(4):
# 		if i==j: continue
# 		print i+j

# print index_mapper(3,1)
# print index_mapper(1,3)

# exit(1)

header =  "nu    "
header += "mu    "
header += "rho   "
header += "sigma "
header += "map(mu, nu):     "
header += "map(rho, sigma): "
header += "lambda-sum"
print header

husk å invertere*(-1) hver gang orginal rekkefølge på indekser er switcha

mu = 0
for nu in range(1,4):
	lambda_indexes = []
	for iLambda in range(0,4):
		if iLambda != mu and iLambda != nu:
			lambda_indexes.append("%d %d %d | %d %d %d" % (mu, iLambda, index_mapper(mu, iLambda), nu, iLambda, index_mapper(nu, iLambda)))

	rho = next_index(nu)
	sigma = next_index(rho)


	# msg  = "%-6d" % mu
	# msg += "%-6d" % nu
	# msg += "%-6d" % rho
	# msg += "%-6d" % sigma
	# msg += "%-17d" % index_mapper(mu, nu)
	# msg += "%-17d" % index_mapper(rho, sigma)
	# msg += " ".join(lambda_indexes)
	# print msg

	print mu, nu, rho, sigma, "map1: %d" % index_mapper(mu, nu), "map2: %d" % index_mapper(rho, sigma), "lambda_indexes: ", lambda_indexes


