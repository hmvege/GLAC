import numpy as np
import matplotlib.pyplot as plt

convert_to_array = lambda s: np.asarray([float(i) for i in s.split(" ")])

table = [
	convert_to_array("0.701 1.75 0.273 0.39"),
	convert_to_array("1.617 8.15 0.30 0.187"),
	convert_to_array("2.244 15.50 0.40 0.177"),
	convert_to_array("3.028 28.14 0.63 0.209"),
	convert_to_array("3.982 48.38 0.81 0.202"),
	convert_to_array("5.167 80.90 0.81 0.157"),
	convert_to_array("1.699 9.07 0.41 0.24"),
	convert_to_array("3.686 41.6 0.83 0.22"),
	convert_to_array("1.750 9.58 0.39 0.22"),
	convert_to_array("3.523 37.8 0.56 0.16"),
	convert_to_array("1.741 9.44 0.35 0.20"),
	convert_to_array("3.266 32.7 0.68 0.21"),
]

def print_values(lattice):
	Q2, Q4, Q4C, R = lattice
	print ("Q2: %10.4f Q4: %10.4f Q4C: %10.4f R: %10.4f | Q4C: %10.4f R: %10.4f" 
		% (Q2, Q4, Q4C, R, Q4 - 3*Q2**2, (Q4 - 3*Q2**2)/Q2))

for lattice in table:
	print_values(lattice)
