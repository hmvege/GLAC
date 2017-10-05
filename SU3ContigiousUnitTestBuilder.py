import numpy as np

U1 = np.asarray([[   1 +    1j ,     1 +    2j ,     1 +    3j   ],
[   2 +    1j ,     2 +    2j ,     2 +    3j   ],
[   3 +    1j ,     3 +    2j ,     3 +    3j   ]],dtype=complex)

U2 = np.asarray([[   4 +    4j ,     4 +    5j ,     4 +    6j  ],
[   5 +    4j ,     5 +    5j ,     5 +    6j  ],
[   6 +    4j ,     6 +    5j ,     6 +    6j  ]],dtype=complex)


print U1,'\n'

print U2,'\n'

print U1*U2