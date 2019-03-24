import numpy as np
import sys
import warnings
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	import matplotlib.pyplot as plt


def get_lattice_spacing(beta):
	"""
	Function for getting the lattice spacing. From paper by Guagnelli et. al., 
	Precision computation of a low-energy reference scale in quenched lattice
	LQCD, 1998.

	Args:
		beta: beta value

	Returns:
		a: lattice spacing in fermi
	"""

	r0 = 0.5
	beta_low = 5.7
	beta_high = 6.57

	if np.any(beta < beta_low) or np.any(beta > beta_high):
		raise Warning("Beta value of %f is outside of defined area [%f, %f]."
			% (beta, beta_low, beta_high))

	def _get_a(b):
		"""Gets the beta value without any error."""
		bval = (b - 6.0)
		_a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*r0
		return _a

	a = _get_a(beta)
	a_err_slope = ((0.6 - 0.3)/100.0)/(beta_high - beta_low) # err% / beta
	a_err_const = 0.3/100 - a_err_slope*beta_low
	a_err_percent = lambda _b: a_err_slope*_b + a_err_const
	a_err = a*a_err_percent(beta)

	return a, a*a_err_percent(beta) # fermi

def main(args):
	if len(args) != 0:
		beta = float(args[0])
		a, a_err = get_lattice_spacing(beta)
		print "beta        = %.2f" % beta
		print "a           = %.8f fm" % a
		print "a_err       = %.8f fm" % a_err
		if len(args) == 2:
			N = float(args[1])
			print "N           = %d" % N
			print "Side length = %.3f +/- %.5f fm" % (N*a, N*a_err)
			print "Volume      = %.3sf +/- %.5f fm^3" % ((N*a)**3, 3*N*a_err*(N*a)**2)
	else:
		beta = 6.0
		a, a_err = get_lattice_spacing(beta)
		N = 24
		print "beta        = %.2f" % beta
		print "a           = %.8f fm" % a
		print "a_err       = %.8f fm" % a_err
		print "N           = %d" % N
		print "Side length = %.3f +/- %.5f fm" % (N*a, N*a_err)
		print "Volume      = %.3f +/- %.5f fm^3" % ((N*a)**3, 3*N*a_err*(N*a)**2)


	if args[-1] == "p":
		X = np.linspace(5.7,6.5,1000)
		Y, Y_err = get_lattice_spacing(X)

		plt.plot(X,Y,label=r"$a = \frac{1}{2}\exp(-1.6805 - 1.7139\beta + 0.8155\beta^2 - 0.6667\beta^3)$")
		plt.axvline(beta,alpha=0.3,color="0",label=r"$\beta=%.2f$" % beta)
		plt.xlabel(r"$\beta$")
		plt.ylabel(r"$a[fm]$")
		plt.title("Beta value vs lattice spacing")
		plt.legend()
		plt.show()
	else:
		print "\nAdd p in args to plot."

if __name__ == '__main__':
	main(sys.argv[1:])