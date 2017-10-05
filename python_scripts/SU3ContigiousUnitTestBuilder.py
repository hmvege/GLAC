import numpy as np

U1 = [
[1 + 1j,              1 + 2j,              1 + 3j],
[2 + 1j,              2 + 2j,              2 + 3j],
[3 + 1j,              3 + 2j,              3 + 3j]]
U1 = np.asarray(U1,dtype=complex)

U2 = [
[4 + 4j,              4 + 5j,              4 + 6j],
[5 + 4j,              5 + 5j,              5 + 6j],
[6 + 4j,              6 + 5j,              6 + 6j]]
U2 = np.asarray(U2,dtype=complex)


print U1,'\n'

print U2,'\n'

print np.dot(U1,U2), '\n'
