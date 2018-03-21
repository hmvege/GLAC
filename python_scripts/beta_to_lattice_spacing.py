import numpy as np, matplotlib.pyplot as plt, sys

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a

def main(args):
	if len(args) != 0:
		beta = float(args[0])
		print "beta        = %.2f" % beta
		print "a           = %.8f fm" % getLatticeSpacing(beta)
		if len(args) == 2:
			N = float(args[1])
			print "N           = %d" % N
			print "Side length = %.2f fm" % (N*getLatticeSpacing(beta))
			print "Volume      = %.2f fm^3" % (N*getLatticeSpacing(beta))**3
	else:
		beta = 6.0
		N = 24
		print "beta        = %.2f" % beta
		print "a           = %.8f fm" % getLatticeSpacing(beta)
		print "N           = %d" % N
		print "Side length = %.2f fm" % (N*getLatticeSpacing(beta))
		print "Volume      = %.2f fm^3" % (N*getLatticeSpacing(beta))**3

	X = np.linspace(5.7,6.5,1000)
	Y = getLatticeSpacing(X)

	plt.plot(X,Y,label=r"$a = \frac{1}{2}\exp(-1.6805 - 1.7139\beta + 0.8155\beta^2 - 0.6667\beta^3)$")
	plt.axvline(beta,alpha=0.3,color="0",label=r"$\beta=%.2f$" % beta)
	plt.xlabel(r"$\beta$")
	plt.ylabel(r"$a[fm]$")
	plt.title("Beta value vs lattice spacing")
	plt.legend()
	plt.show()

if __name__ == '__main__':
	main(sys.argv[1:])